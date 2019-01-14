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

!> summary: Software subsystem in charge of implementing the Jacobi 
!> (a.k.a. diagonal scaling) preconditioner operator.

module jacobi_preconditioner_names
 ! Tools
 use types_names
 use list_types_names
 use FPL
 
 ! Integration related modules
 use triangulation_names
 use fe_space_names
 use fe_operator_names
 
 ! Linear Algebra related modules
 use operator_names
 use par_scalar_array_names
 use serial_scalar_array_names
 use vector_space_names
 use vector_names
 
 use matrix_names
 use sparse_matrix_parameters_names
 use sparse_matrix_names
 use par_sparse_matrix_names
 
 ! Parallel communication-related data structures
 use environment_names
 
 implicit none
# include "debug.i90"
 private
 
 integer(ip), parameter :: jacobi_preconditioner_STATE_START    = 0 !< Entity has no stored data
 integer(ip), parameter :: jacobi_preconditioner_STATE_CREATED  = 1 !< Ready to build preconditioner
 integer(ip), parameter :: jacobi_preconditioner_STATE_SYMBOLIC = 2 !< Symbolic data already computed
 integer(ip), parameter :: jacobi_preconditioner_STATE_NUMERIC  = 3 !< Numerical data already computed
 
  !> 
  !> This module contains the software subsystem in charge of implementing the Jacobi 
  !> (a.k.a. diagonal scaling) preconditioner operator. If one aims to solve the linear system 
  !> A x = f at hand of the simplest iterative method, i.e. Richardsond, a Jacobi-preconditioned 
  !> iteration looks like x^{k+1} = x^k + inv(D) * ( A x^k - f ). Here, M = D = diag(A) is the 
  !> preconditioner.
  !>
  !> # Basic usage at the driver level
  !>
  !> ```fortran
  !> ...
  !>    type(jacobi_preconditioner_t) :: jp ! Declaration in the preamble
  !> ...
  !>    ! Setup preconditioner after assembly of linear system
  !>    call jp%create(fe_affine_operator) ! Or jp%create(fe_affine_operator,parameter_list)
  !>    call jp%symbolic_setup()
  !>    call jp%numerical_setup()
  !> ...
  !>    call iterative_linear_solver%set_operators(fe_affine_operator%get_tangent(),jp)
  !> ...
  !>    call jp%free()
  !> ...
  !> ```
  !>
  !> @note To save memory usage, the inverse of the diagonal matrix of the FE matrix, i.e. 
  !> inv(D) is actually not stored. As the only nontrivial entries of inv(D) are, precisely, 
  !> the diagonal ones, only these are computed and stored in the [[jacobi_preconditioner_t:inverse_diagonal]]
  !> polymorphic [[vector_t]] instance with size equal to the local FE matrix size (n). 
  !> In other words, instead of inv(D), we actually compute and store {inv(D)_ii}, i = 1,...,n. 
  !> Matrix-vector product is then performed as an entrywise product, i.e. y = inv(D) * x is 
  !> computed entry by entry as y_i = inv(D)_{ii} * x_i, i = 1,...,n.
  !> 
  !> @warning [[jacobi_preconditioner_t]] assumes the FE matrix is positive definite 
  !> and has a single-block structure. Execution stops at run-time whenever any of 
  !> these conditions is not verified.
  !>
  !> # State transition diagram for type [[jacobi_preconditioner_t]]
  !> 
  !> | Input state | Action | Output state | Comments |
  !> | ----------- | :----: | :----------: | :------: |
  !> | Start       | create                                | Created  | |
  !> | Start       | free_clean                            | Start    | |
  !> | Start       | free_symbolic                         | Start    | |
  !> | Start       | free_numeric                          | Start    | |
  !> | Start       | update_matrix                         | Start    | it does nothing |
  !> | Created     | symbolic_setup                        | Symbolic | symbolic_setup() |
  !> | Created     | numerical_setup                       | Numeric  | symbolic_setup()+numerical_setup() |
  !> | Created     | apply                                 | Numeric  | symbolic_setup()+numerical_setup() |
  !> | Created     | free_clean                            | Start    | |
  !> | Created     | free_symbolic                         | Create   | it does nothing |
  !> | Created     | free_numeric                          | Create   | it does nothing |
  !> | Created     | update_matrix                         | Create   | it does nothing |
  !> | Symbolic    | symbolic_setup                        | Symbolic | it does nothing |
  !> | Symbolic    | numerical_setup                       | Numeric  | numerical_setup() |
  !> | Symbolic    | apply                                 | Numeric  | numerical_setup() |
  !> | Symbolic    | free_clean                            | Start    | |
  !> | Symbolic    | free_symbolic                         | Created  | |
  !> | Symbolic    | free_numeric                          | Symbolic | it does nothing |
  !> | Symbolic    | update_matrix + same_nonzero_pattern  | Symbolic | it does nothing |
  !> | Symbolic    | update_matrix + !same_nonzero_pattern | Symbolic | free_symbolic()+symbolic_setup() |
  !> | Numeric     | symbolic_setup                        | Numeric  | it does nothing |
  !> | Numeric     | numeric_setup                         | Numeric  | it does nothing |
  !> | Numeric     | apply                                 | Numeric  | it does nothing |
  !> | Numeric     | free_numeric                          | Symbolic | |
  !> | Numeric     | free_symbolic                         | Created  | |
  !> | Numeric     | free_clean                            | Start    | |
  !> | Numeric     | update_matrix + same_nonzero_pattern  | Numeric  | free_numerical_setup()+numerical_setup() |
  !> | Numeric     | update_matrix + !same_nonzero_pattern | Numeric  | free_numerical_setup()+free_symbolic_setup()+symbolic_setup()+numeric_setup()
  !>
type, extends(operator_t) :: jacobi_preconditioner_t
   private
   integer(ip)                        :: state = jacobi_preconditioner_STATE_START !< Refer to the state diagramme of [[jacobi_preconditioner_t]]
   class(environment_t) , pointer     :: environment           => NULL()
   type(fe_operator_t)  , pointer     :: fe_nonlinear_operator => NULL()           !< Pointer to the [[fe_operator_t]] this [[jacobi_preconditioner_t]] instance has been created from.
   class(vector_t)      , allocatable :: inverse_diagonal        
   type(parameterlist_t), pointer     :: jacobi_preconditioner_params => NULL()    !< Pointer to a [[parameter_list_t]] to customize the preconditioner (e.g. set a relaxation parameter).
   !< @warning The [[jacobi_preconditioner_t:jacobi_preconditioner_params]] pointer is set-up 
   !< during [[jacobi_preconditioner_t:create]] and re-used in the rest of stages. Therefore,
   !< [[parameter_list_t]] to which [[jacobi_preconditioner_t]] points to MUST NOT BE freed  
   !< before [[jacobi_preconditioner_t]].
 contains
 
   ! Creational methods
   procedure, non_overridable  :: jacobi_preconditioner_create_w_parameter_list
   procedure, non_overridable  :: jacobi_preconditioner_create_wo_parameter_list
   generic                     :: create => jacobi_preconditioner_create_w_parameter_list, &
                                            jacobi_preconditioner_create_wo_parameter_list
   
   procedure,                  private :: create_vector_spaces   => jacobi_preconditioner_create_vector_spaces
   procedure,                  private :: create_and_allocate_inverse_diagonal & 
                                                                 => jacobi_preconditioner_create_and_allocate_inverse_diagonal
 
   ! State transition handling-related TBPs
   procedure, non_overridable, private :: set_state_start        => jacobi_preconditioner_set_state_start
   procedure, non_overridable, private :: set_state_created      => jacobi_preconditioner_set_state_created
   procedure, non_overridable, private :: set_state_symbolic     => jacobi_preconditioner_set_state_symbolic
   procedure, non_overridable, private :: set_state_numeric      => jacobi_preconditioner_set_state_numeric
   procedure, non_overridable, private :: state_is_start         => jacobi_preconditioner_state_is_start
   procedure, non_overridable, private :: state_is_created       => jacobi_preconditioner_state_is_created
   procedure, non_overridable, private :: state_is_symbolic      => jacobi_preconditioner_state_is_symbolic
   procedure, non_overridable, private :: state_is_numeric       => jacobi_preconditioner_state_is_numeric
 
   ! Symbolic setup-related TBPs
   procedure, non_overridable          :: symbolic_setup         => jacobi_preconditioner_symbolic_setup

   ! Numerical setup-related TBPs
   procedure, non_overridable          :: numerical_setup        => jacobi_preconditioner_numerical_setup

   ! Apply related TBPs
   procedure                           :: apply                  => jacobi_preconditioner_apply
   procedure                           :: apply_add              => jacobi_preconditioner_apply_add
   
   ! Free-related TBPs
   procedure, non_overridable          :: free                   => jacobi_preconditioner_free
   procedure, non_overridable          :: free_clean             => jacobi_preconditioner_free_clean
   procedure, non_overridable          :: free_symbolic_setup    => jacobi_preconditioner_free_symbolic_setup   
   procedure, non_overridable          :: free_numerical_setup   => jacobi_preconditioner_free_numerical_setup
   procedure, non_overridable, private :: free_and_destroy_inverse_diagonal & 
                                                                 => jacobi_preconditioner_free_and_destroy_inverse_diagonal
 
   procedure, non_overridable, private :: am_i_l1_task           => jacobi_preconditioner_am_i_l1_task
   procedure                           :: is_linear              => jacobi_preconditioner_is_linear
   procedure, private                  :: get_par_environment    => jacobi_preconditioner_get_par_environment
   procedure, private                  :: set_par_environment    => jacobi_preconditioner_set_par_environment
   
   ! Miscellaneous 
   procedure, private                  :: get_par_sparse_matrix  => jacobi_preconditioner_get_par_sparse_matrix
   procedure, private                  :: get_fe_space           => jacobi_preconditioner_get_fe_space
   procedure, private                  :: get_par_fe_space       => jacobi_preconditioner_get_par_fe_space
   procedure                 , private :: is_operator_associated => jacobi_preconditioner_is_operator_associated
   procedure                 , private :: nullify_operator       => jacobi_preconditioner_nullify_operator 
     
   ! Update-matrix related TBPs
   procedure                           :: update_matrix          => jacobi_preconditioner_update_matrix
end type jacobi_preconditioner_t

 ! See issue #270 (open as of Jan 11th 2019)
 ! Cannot use constructor s.t. instead of following the basic usage, i.e. jp%create(), 
 ! jp%symbolic_setup(), ..., ils%set_operators(fe_op,jp), we straighforwardly set the 
 ! operators with a constructor, i.e. ils%set_operators(fe_op,jp_t(fe_op)), because the 
 ! jp instance must be temporal and assigned to an [[lvalue_operator_t]]. At the current 
 ! design of [[lvalue_operator_t]] (commit ID 11170497), both requirements are incompatible 
 ! because [[lvalue_operator_t]] can only assign a permanent polymorphic instance of 
 ! [[fe_operator_t]].
 !interface jacobi_preconditioner_t
 !  module procedure create_jacobi_preconditioner
 !end interface jacobi_preconditioner_t
 
 public :: jacobi_preconditioner_t

contains

 subroutine jacobi_preconditioner_create_w_parameter_list ( this, fe_nonlinear_operator, jacobi_preconditioner_params )
   implicit none
   class(jacobi_preconditioner_t)        , intent(inout) :: this
   class(fe_operator_t)          , target, intent(in)    :: fe_nonlinear_operator
   type(parameterlist_t)         , target, intent(in)    :: jacobi_preconditioner_params
   call this%create(fe_nonlinear_operator)
   this%jacobi_preconditioner_params => jacobi_preconditioner_params
 end subroutine jacobi_preconditioner_create_w_parameter_list

 subroutine jacobi_preconditioner_create_wo_parameter_list ( this, fe_nonlinear_operator )
   implicit none
   class(jacobi_preconditioner_t)        , intent(inout) :: this
   class(fe_operator_t)          , target, intent(in)    :: fe_nonlinear_operator
   class(matrix_t)       , pointer :: fe_matrix
   type(par_fe_space_t)  , pointer :: fe_space
   class(triangulation_t), pointer :: triangulation
   call this%free()
   assert ( this%state_is_start() )
   this%fe_nonlinear_operator => fe_nonlinear_operator
   fe_space => this%get_par_fe_space()
   triangulation => fe_space%get_triangulation()
   massert ( fe_space%get_num_blocks() == 1, 'jacobi_preconditioner_create: FE matrix must have a single block' )
   call this%set_par_environment(triangulation%get_environment())
   if ( this%am_i_l1_task() ) then 
     fe_matrix => this%fe_nonlinear_operator%get_matrix()
     mcheck( fe_matrix%get_sign() == SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE, 'jacobi_preconditioner_create: FE matrix must be positive definite' )
   end if
   call this%create_vector_spaces()
   call this%set_state_created()
   nullify(this%jacobi_preconditioner_params)
 end subroutine jacobi_preconditioner_create_wo_parameter_list

 subroutine jacobi_preconditioner_create_vector_spaces (this)
   implicit none
   class(jacobi_preconditioner_t), intent(inout)  :: this
   type(vector_space_t), pointer :: fe_nonlinear_operator_domain_vector_space
   type(vector_space_t), pointer :: fe_nonlinear_operator_range_vector_space
   type(vector_space_t), pointer :: jacobi_preconditioner_domain_vector_space
   type(vector_space_t), pointer :: jacobi_preconditioner_range_vector_space
   fe_nonlinear_operator_domain_vector_space => this%fe_nonlinear_operator%get_domain_vector_space()
   fe_nonlinear_operator_range_vector_space  => this%fe_nonlinear_operator%get_range_vector_space()
   assert ( fe_nonlinear_operator_domain_vector_space%equal_to(fe_nonlinear_operator_range_vector_space) )
   jacobi_preconditioner_domain_vector_space => this%get_domain_vector_space()
   jacobi_preconditioner_range_vector_space  => this%get_range_vector_space()
   call fe_nonlinear_operator_domain_vector_space%clone(jacobi_preconditioner_domain_vector_space)
   call fe_nonlinear_operator_range_vector_space%clone(jacobi_preconditioner_range_vector_space)
   !< @note The type [[fe_nonlineal_operator_t]] must have equal domain and   
   !< range vector spaces. In other words, the FE matrix must be squared.
 end subroutine jacobi_preconditioner_create_vector_spaces

 !> summary: Creates [[jacobi_preconditioner:inverse_diagonal]]
 !> with memory space for the diagonal entries of inv(D) .
 subroutine jacobi_preconditioner_create_and_allocate_inverse_diagonal ( this ) 
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   assert ( this%state_is_created() .or. this%state_is_symbolic() .or. this%state_is_numeric() )
   call this%free_and_destroy_inverse_diagonal()
   call this%create_domain_vector(this%inverse_diagonal)
   !< [[jacobi_preconditioner_t:inverse_diagonal]] is a polymorphic [[vector_t]] 
   !< instance with size equal to the local FE matrix size (n). To save memory usage,
   !< [[jacobi_preconditioner_t:inverse_diagonal]] only stores the diagonal entries 
   !< of inv(D). Since D is a diagonal matrix (D_{ij} = 0, i != j), the matrix-vector product   
   !< inv(D)*x, with x a given vector in the FE nonlinear operator domain=range space, 
   !< can be reformulated as an entrywise vector product between the diagonal entries of 
   !< inv(D) and x, i.e. y = inv(D)*x satisfies y_i = inv(D)_{ii} * x_i, i = 1,...,n.
 end subroutine jacobi_preconditioner_create_and_allocate_inverse_diagonal 

 subroutine jacobi_preconditioner_set_state_start(this)
   class(jacobi_preconditioner_t), intent(inout) :: this
   this%state = jacobi_preconditioner_STATE_START
 end subroutine jacobi_preconditioner_set_state_start

 subroutine jacobi_preconditioner_set_state_created(this)
   class(jacobi_preconditioner_t), intent(inout) :: this
   this%state = jacobi_preconditioner_STATE_CREATED
 end subroutine jacobi_preconditioner_set_state_created

 subroutine jacobi_preconditioner_set_state_symbolic(this)
   class(jacobi_preconditioner_t), intent(inout) :: this
   this%state = jacobi_preconditioner_STATE_SYMBOLIC
 end subroutine jacobi_preconditioner_set_state_symbolic

 subroutine jacobi_preconditioner_set_state_numeric(this)
   class(jacobi_preconditioner_t), intent(inout) :: this
   this%state = jacobi_preconditioner_STATE_NUMERIC
 end subroutine jacobi_preconditioner_set_state_numeric

 function jacobi_preconditioner_state_is_start(this) result(is_start)
   class(jacobi_preconditioner_t), intent(in) :: this
   logical                                 :: is_start
   is_start = this%state == jacobi_preconditioner_STATE_START
 end function jacobi_preconditioner_state_is_start

 function jacobi_preconditioner_state_is_created(this) result(is_start)
   class(jacobi_preconditioner_t), intent(in) :: this
   logical                                 :: is_start
   is_start = this%state == jacobi_preconditioner_STATE_CREATED
 end function jacobi_preconditioner_state_is_created

 function jacobi_preconditioner_state_is_symbolic(this) result(is_symbolic_setup)
   class(jacobi_preconditioner_t), intent(in) :: this
   logical                                 :: is_symbolic_setup
   is_symbolic_setup = this%state == jacobi_preconditioner_STATE_SYMBOLIC
 end function jacobi_preconditioner_state_is_symbolic

 function jacobi_preconditioner_state_is_numeric(this) result(is_numerical_setup)
   class(jacobi_preconditioner_t), intent(in) :: this
   logical                                 :: is_numerical_setup
   is_numerical_setup= this%state == jacobi_preconditioner_STATE_NUMERIC
 end function jacobi_preconditioner_state_is_numeric

 !> summary: Setup phase in charge of creating and allocating [[jacobi_preconditioner:inverse_diagonal]].
 !> See also its main called TBP [[jacobi_preconditioner:create_and_allocate_inverse_diagonal]].
 subroutine jacobi_preconditioner_symbolic_setup ( this )
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   type(environment_t), pointer :: par_environment  
   par_environment => this%get_par_environment()
   assert ( this%state_is_created() .or. this%state_is_symbolic() .or. this%state_is_numeric() )
   if ( this%state_is_created() ) then
     call this%create_vector_spaces()
     !< In the case the call to this subroutine is triggered by 
     !< [[jacobi_preconditioner_t:update_matrix]] with 
     !< 'same_nonzero_pattern=.false.' it might be necessary to 
     !< re-generate vector spaces associated to this. Provided that 
     !< we do not know who triggered this subroutine, and that the 
     !< computational time spent here is not significant, we always 
     !< regenerate vector spaces right before returning control
     !< from a call to symbolic_setup.
     call this%create_and_allocate_inverse_diagonal()
     call this%set_state_symbolic()
   end if
 end subroutine jacobi_preconditioner_symbolic_setup

 !> summary: Setup phase in charge of filling the *globally-assembled*
 !> entries of [[jacobi_preconditioner:inverse_diagonal]].
 subroutine jacobi_preconditioner_numerical_setup ( this )
   implicit none
   class(jacobi_preconditioner_t), intent(inout)   :: this
   type(environment_t)        , pointer :: par_environment
   class(matrix_t)            , pointer :: fe_matrix
   type(serial_scalar_array_t), pointer :: serial_scalar_array
   real(rp)                   , pointer :: id_entries(:)
   assert ( this%state_is_created() .or. this%state_is_symbolic() .or. this%state_is_numeric() )
   if ( this%state_is_created() ) then
     call this%symbolic_setup()
   end if 
   assert ( this%state_is_symbolic() .or. this%state_is_numeric() )
   if ( this%state_is_symbolic() ) then
     par_environment => this%get_par_environment()
     if ( par_environment%am_i_l1_task() ) then
       fe_matrix => this%fe_nonlinear_operator%get_matrix()
       select type ( id_vector => this%inverse_diagonal )
         class is (serial_scalar_array_t)
           id_entries          => id_vector%get_entries()
           call fe_matrix%extract_diagonal(id_entries)
         class is (par_scalar_array_t)
           serial_scalar_array => id_vector%get_serial_scalar_array()
           id_entries          => serial_scalar_array%get_entries()
           call fe_matrix%extract_diagonal(id_entries)
         class default
           check(.false.)
       end select
     end if
     call this%inverse_diagonal%comm()
     call this%inverse_diagonal%entrywise_invert()
     call this%set_state_numeric()
   end if
 end subroutine jacobi_preconditioner_numerical_setup

 !> summary: Apply matrix-vector product y = inv(D)*x
 subroutine jacobi_preconditioner_apply (this, x, y)
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   class(vector_t)               , intent(in)    :: x
   class(vector_t)               , intent(inout) :: y
   assert ( this%state_is_created() .or. this%state_is_symbolic() .or. this%state_is_numeric() )
   if (this%state_is_created() ) then
     call this%symbolic_setup()
   end if 
   assert ( this%state_is_symbolic() .or. this%state_is_numeric() )
   if (this%state_is_symbolic() ) then
     call this%numerical_setup()
   end if 
   assert ( this%state_is_numeric() )
   call this%abort_if_not_in_domain(x)
   call this%abort_if_not_in_range(y)
   call x%GuardTemp()
   call y%entrywise_product(this%inverse_diagonal,x)
   !< y = inv(D)*x satisfies y_i = inv(D)_{ii} * x_i, i = 1,...,n, because D_{ij} = 0, i != j.
   !< See also [[jacobi_preconditioner:create_and_allocate_inverse_diagonal]].
   call x%CleanTemp()
 end subroutine jacobi_preconditioner_apply

 !> summary: Apply matrix-vector product y = inv(D)*x + y
 subroutine jacobi_preconditioner_apply_add(this, x, y)
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   class(vector_t)               , intent(in)    :: x
   class(vector_t)               , intent(inout) :: y
   class(vector_t)     , allocatable :: w
   type(vector_space_t), pointer     :: range_vector_space
   integer(ip)                       :: istat
   call this%abort_if_not_in_domain(x)
   call this%abort_if_not_in_range(y)
   call x%GuardTemp()
   range_vector_space => this%get_range_vector_space()
   call range_vector_space%create_vector(w)
   call this%apply(x,w)
   call y%axpby(1.0, w, 1.0)
   call x%CleanTemp()
   call w%free()
   deallocate(w, stat=istat); check(istat==0)
 end subroutine jacobi_preconditioner_apply_add

 subroutine jacobi_preconditioner_free(this)
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   call this%free_numerical_setup()
   call this%free_symbolic_setup()
   call this%free_clean()
 end subroutine jacobi_preconditioner_free

 !> summary: destroys all memory allocated by or pointers associated by the member 
 !> variables of [[jacobi_preconditioner_t]]
 subroutine jacobi_preconditioner_free_clean(this)
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   assert ( this%state_is_start() .or. this%state_is_created() .or. this%state_is_symbolic() .or. this%state_is_numeric() )
   if ( this%state_is_numeric() ) then
     call this%free_numerical_setup() 
   end if 
   if ( this%state_is_symbolic() ) then
     call this%free_symbolic_setup() 
   end if
   call this%nullify_operator()
   call this%free_vector_spaces()
   nullify(this%jacobi_preconditioner_params)
   nullify(this%environment)
   call this%set_state_start()
 end subroutine jacobi_preconditioner_free_clean

 subroutine jacobi_preconditioner_free_symbolic_setup(this)
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   assert ( this%state_is_start() .or. this%state_is_created() .or. this%state_is_symbolic() .or. this%state_is_numeric() )
   if ( this%state_is_numeric() ) then
     call this%free_numerical_setup() 
   end if
   if ( this%state_is_symbolic() ) then 
     call this%free_and_destroy_inverse_diagonal()
     call this%set_state_created()
   end if
 end subroutine jacobi_preconditioner_free_symbolic_setup

 subroutine jacobi_preconditioner_free_numerical_setup(this)
   implicit none
   class(jacobi_preconditioner_t)           , intent(inout) :: this
   type(environment_t)   , pointer       :: par_environment
   assert ( this%state_is_start() .or. this%state_is_created() .or. this%state_is_symbolic() .or. this%state_is_numeric() )
   if ( this%state_is_numeric() ) then
     call this%inverse_diagonal%init(0.0_rp)
     call this%set_state_symbolic()
   end if
 end subroutine jacobi_preconditioner_free_numerical_setup

 subroutine jacobi_preconditioner_free_and_destroy_inverse_diagonal ( this ) 
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   integer(ip) :: istat
   if ( allocated(this%inverse_diagonal) ) then
     call this%inverse_diagonal%free()
     deallocate(this%inverse_diagonal, stat=istat); check(istat==0);
   end if
 end subroutine jacobi_preconditioner_free_and_destroy_inverse_diagonal 

 function jacobi_preconditioner_get_par_environment(this)
   implicit none
   class(jacobi_preconditioner_t), target, intent(in) :: this
   class(environment_t), pointer :: jacobi_preconditioner_get_par_environment
   jacobi_preconditioner_get_par_environment => this%environment
 end function jacobi_preconditioner_get_par_environment

 subroutine jacobi_preconditioner_set_par_environment(this, environment)
   implicit none
   class(jacobi_preconditioner_t)        , intent(inout) :: this
   class(environment_t)          , target, intent(in)    :: environment
   this%environment => environment
 end subroutine jacobi_preconditioner_set_par_environment

 function jacobi_preconditioner_get_par_sparse_matrix(this)
   !< Helper function that extracts a run-time polymorphic [[matrix_t]]
   !< from the [[fe_operator_t]], and dynamically casts it into  
   !< type [[par_sparse_matrix_t]]. If the dynamic cast cannot be performed, 
   !< because class [[matrix_t]] is NOT of type [[par_sparse_matrix_t], then  
   !< it aborts the execution of the program.
   implicit none
   class(jacobi_preconditioner_t), intent(in) :: this
   type(par_sparse_matrix_t), pointer :: jacobi_preconditioner_get_par_sparse_matrix
   class(matrix_t)          , pointer :: matrix
   matrix => this%fe_nonlinear_operator%get_matrix()
   select type( matrix )
   type is ( par_sparse_matrix_t )
      jacobi_preconditioner_get_par_sparse_matrix => matrix
      class default
      check(.false.)
   end select
 end function jacobi_preconditioner_get_par_sparse_matrix

 function jacobi_preconditioner_get_fe_space(this)
   implicit none
   class(jacobi_preconditioner_t), intent(in) :: this
   class(base_fe_space_t), pointer :: jacobi_preconditioner_get_fe_space
   jacobi_preconditioner_get_fe_space => this%get_par_fe_space()
 end function jacobi_preconditioner_get_fe_space

 function jacobi_preconditioner_get_par_fe_space(this)
   !< Helper function that extracts a run-time polymorphic [[serial_fe_space_t]]
   !< from the [[fe_operator_t]], and dynamically casts it into  
   !< type [[par_fe_space_t]]. If the dynamic cast cannot be performed,
   !< because class [[serial_fe_space_t]] is NOT of type [[par_fe_space_t]],  
   !< then it aborts the execution of the program.
   implicit none
   class(jacobi_preconditioner_t), intent(in) :: this
   type(par_fe_space_t)    , pointer    :: jacobi_preconditioner_get_par_fe_space
   class(serial_fe_space_t), pointer    :: fe_space
   fe_space => this%fe_nonlinear_operator%get_fe_space()
   select type(fe_space)
   class is (par_fe_space_t)
      jacobi_preconditioner_get_par_fe_space => fe_space
      class default
      check(.false.)
   end select
 end function jacobi_preconditioner_get_par_fe_space

 function jacobi_preconditioner_am_i_l1_task(this)
   implicit none
   class(jacobi_preconditioner_t), intent(in) :: this
   logical :: jacobi_preconditioner_am_i_l1_task
   type(environment_t), pointer     :: par_environment
   par_environment => this%get_par_environment()
   jacobi_preconditioner_am_i_l1_task = par_environment%am_i_l1_task()
 end function jacobi_preconditioner_am_i_l1_task

 function jacobi_preconditioner_is_linear( this )
   implicit none
   class(jacobi_preconditioner_t), intent(in) :: this
   logical :: jacobi_preconditioner_is_linear
   jacobi_preconditioner_is_linear = .true.
 end function jacobi_preconditioner_is_linear

 function jacobi_preconditioner_is_operator_associated( this )
   implicit none
   class(jacobi_preconditioner_t), intent(in) :: this
   logical :: jacobi_preconditioner_is_operator_associated
   jacobi_preconditioner_is_operator_associated = associated(this%fe_nonlinear_operator)
 end function jacobi_preconditioner_is_operator_associated

 subroutine jacobi_preconditioner_nullify_operator ( this )
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   nullify(this%fe_nonlinear_operator)
 end subroutine jacobi_preconditioner_nullify_operator 

 !> summary: Updates the [[jacobi_preconditioner_t]] instance, specifically, the 
 !> entries of the inverted diagonal of the FE matrix, whenever FE matrix has been updated.
 subroutine jacobi_preconditioner_update_matrix(this, same_nonzero_pattern )
   !< This TBP should be (generally) called right after the FE matrix has been updated.
   !< It is in charge of updating the [[jacobi_preconditioner_t:inverse_diagonal]] s.t.
   !< it matches the new modified FE matrix. If the FE matrix preserves its nonzero pattern,
   !< numerical setup is performed. Ow, symbolic setup is also performed. Actually, symbolic 
   !< setup is only necessary when the num rows (= num columns) changes, but this optimization
   !< would probably not lead to significant performance improvement.
   implicit none
   class(jacobi_preconditioner_t), intent(inout) :: this
   logical                       , intent(in)    :: same_nonzero_pattern
   assert ( this%state_is_created() .or. this%state_is_symbolic() .or. this%state_is_numeric() )
   if ( same_nonzero_pattern ) then
     if ( this%state_is_numeric() ) then
       call this%free_numerical_setup()
       call this%numerical_setup()
     end if   
   else
     if ( this%state_is_numeric() ) then
       call this%free_numerical_setup()
       call this%free_symbolic_setup()
       call this%symbolic_setup()
       call this%numerical_setup()
     else if ( this%state_is_symbolic() ) then
       call this%free_symbolic_setup()
       call this%symbolic_setup()
     end if
   end if 
 end subroutine jacobi_preconditioner_update_matrix

 !> summary: Creates a *temporary* instance of [[jacobi_preconditioner_t]]
 !> from an instance of [[fe_operator_t]].
 function create_jacobi_preconditioner(op) result (res)
   implicit none
   class(fe_operator_t)         , intent(in)  :: op
   type(jacobi_preconditioner_t)              :: res
   call op%GuardTemp()
   call res%create(op)
   call res%setTemp()
   call op%CleanTemp()
 end function create_jacobi_preconditioner

end module jacobi_preconditioner_names
