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
module par_preconditioner_dd_diagonal_names
  ! Serial modules
  use types_names
  use memor_names
  use serial_scalar_array_names
  use preconditioner_names, only: invert_diagonal, apply_diagonal, extract_diagonal

  ! Parallel modules
  use par_scalar_array_names
  use par_scalar_matrix_names
  use par_context_names
  use par_environment_names
  use dof_distribution_names
  use psb_penv_mod_names

  ! Abstract modules
  use matrix_names
  use vector_names
  use operator_names
  use fe_affine_operator_names
  
  implicit none
# include "debug.i90"

  private

  type, extends(operator_t) :: par_preconditioner_dd_diagonal_t
     ! Reference to parallel matrix
     type( par_scalar_matrix_t ), pointer     :: p_mat => NULL()   
     real(rp)          , allocatable :: d(:)            ! Inverse of main diagonal
   contains
     procedure :: is_linear => par_preconditioner_dd_diagonal_is_linear
     procedure :: apply     => par_preconditioner_dd_diagonal_apply
     procedure :: free      => par_preconditioner_dd_diagonal_free
  end type par_preconditioner_dd_diagonal_t
  
  interface par_preconditioner_dd_diagonal_create
    module procedure par_preconditioner_dd_diagonal_create_w_par_matrix, par_preconditioner_dd_diagonal_create_w_fe_affine_operator
  end interface
  
  
  ! Types
  public :: par_preconditioner_dd_diagonal_t

  ! Functions
  public :: par_preconditioner_dd_diagonal_create, par_preconditioner_dd_diagonal_free, &
            par_preconditioner_dd_diagonal_symbolic_setup, par_preconditioner_dd_diagonal_numerical_setup, &
            par_preconditioner_dd_diagonal_apply_all_unk

  contains


  !=============================================================================
  subroutine par_preconditioner_dd_diagonal_create_w_par_matrix (p_matrix, p_prec_dd_diagonal)
    implicit none
    ! Parameters
    type(par_scalar_matrix_t)             , target, intent(in)  :: p_matrix
    type(par_preconditioner_dd_diagonal_t)        , intent(out) :: p_prec_dd_diagonal

    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist_domain) )

    p_prec_dd_diagonal%p_mat    => p_matrix
  end subroutine par_preconditioner_dd_diagonal_create_w_par_matrix
  
  !=============================================================================
  subroutine par_preconditioner_dd_diagonal_create_w_fe_affine_operator (fe_affine_operator, p_prec_dd_diagonal)
    implicit none
    ! Parameters
    type(fe_affine_operator_t)            , intent(in)  :: fe_affine_operator
    type(par_preconditioner_dd_diagonal_t), intent(out) :: p_prec_dd_diagonal
    
    ! Local variables
     class(matrix_t)          , pointer :: matrix
     type(par_scalar_matrix_t), pointer :: p_mat
  
     matrix => fe_affine_operator%get_matrix()
     select type(matrix)
     class is(par_scalar_matrix_t)
       p_mat => matrix
     class default
       check(.false.)
     end select 

    call par_preconditioner_dd_diagonal_create_w_par_matrix ( p_mat, p_prec_dd_diagonal )
  end subroutine par_preconditioner_dd_diagonal_create_w_fe_affine_operator

  !=============================================================================
  subroutine par_preconditioner_dd_diagonal_symbolic_setup (p_prec_dd_diagonal)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_diagonal_t)        , intent(inout) :: p_prec_dd_diagonal
  end subroutine par_preconditioner_dd_diagonal_symbolic_setup

  !=============================================================================
  subroutine par_preconditioner_dd_diagonal_numerical_setup ( p_prec_dd_diagonal)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_diagonal_t), target, intent(inout) :: p_prec_dd_diagonal

    ! Locals
    type(par_scalar_matrix_t), pointer :: p_matrix
    integer(ip)       :: neq
    type (par_scalar_array_t) :: p_vec

    p_matrix => p_prec_dd_diagonal%p_mat

    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist_domain) )

    if ( p_prec_dd_diagonal%p_mat%p_env%p_context%iam < 0 ) return

    neq = p_matrix%serial_scalar_matrix%graph%nv
       
    ! Allocate, extract, sum, and invert diagonal 

    ! Allocate + extract
    call memalloc ( neq, p_prec_dd_diagonal%d, __FILE__,__LINE__)
    call extract_diagonal ( p_matrix%serial_scalar_matrix%graph%symmetric_storage, &
                            p_matrix%serial_scalar_matrix%graph%nv, &
                            p_matrix%serial_scalar_matrix%graph%ia, &
                            p_matrix%serial_scalar_matrix%graph%ja, &
                            p_matrix%serial_scalar_matrix%a, &
                            p_matrix%serial_scalar_matrix%graph%nv, &
                            p_prec_dd_diagonal%d )
    
    ! Create a view of p_prec_dd_diagonal%d. This is ugly and
    ! violates principles of OO design, but at least it works
    p_vec%dof_dist      => p_prec_dd_diagonal%p_mat%dof_dist_domain
    p_vec%p_env         => p_prec_dd_diagonal%p_mat%p_env
    call p_vec%serial_scalar_array%create(neq)
    call p_vec%serial_scalar_array%set_view_entries(p_prec_dd_diagonal%d)

    ! Communicate
    call p_vec%comm()
    
    ! Invert diagonal
    call invert_diagonal  ( neq, p_prec_dd_diagonal%d )

  end subroutine par_preconditioner_dd_diagonal_numerical_setup

  !=============================================================================
  subroutine par_preconditioner_dd_diagonal_apply_all_unk (p_prec_dd_diagonal, x, y)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_diagonal_t) , intent(in)    :: p_prec_dd_diagonal
    type(par_scalar_array_t)              , intent(in)    :: x
    type(par_scalar_array_t)              , intent(inout) :: y

    assert ( associated(p_prec_dd_diagonal%p_mat) )
    assert ( associated(p_prec_dd_diagonal%p_mat%p_env) )
    assert ( p_prec_dd_diagonal%p_mat%p_env%created )
    assert ( associated(p_prec_dd_diagonal%p_mat%dof_dist_domain) )

    if(p_prec_dd_diagonal%p_mat%p_env%p_context%iam<0) return

    call apply_diagonal  ( p_prec_dd_diagonal%p_mat%serial_scalar_matrix%graph%nv, & 
                           p_prec_dd_diagonal%d, & 
                           x%serial_scalar_array%b, & 
                           y%serial_scalar_array%b )
  end subroutine par_preconditioner_dd_diagonal_apply_all_unk

  !=============================================================================
  subroutine  par_preconditioner_dd_diagonal_free_in_stages (p_prec_dd_diagonal, mode)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_diagonal_t),  intent(inout) :: p_prec_dd_diagonal
    integer(ip)                  ,  intent(in)    :: mode

    assert ( associated(p_prec_dd_diagonal%p_mat) )
    assert ( associated(p_prec_dd_diagonal%p_mat%p_env) )
    assert ( p_prec_dd_diagonal%p_mat%p_env%created )
    assert ( associated(p_prec_dd_diagonal%p_mat%dof_dist_domain) )

    if (p_prec_dd_diagonal%p_mat%p_env%p_context%iam<0) return
    
    if (mode == free_clean) then
       nullify(p_prec_dd_diagonal%p_mat)
    else if (mode == free_values) then
       call memfree ( p_prec_dd_diagonal%d,__FILE__,__LINE__)
    end if

  end subroutine par_preconditioner_dd_diagonal_free_in_stages

  !=============================================================================
  subroutine par_preconditioner_dd_diagonal_apply (op, x, y)
    implicit none
    ! Parameters
    class(par_preconditioner_dd_diagonal_t)    , intent(in)    :: op
    class(vector_t)   , intent(in)    :: x
    class(vector_t)   , intent(inout) :: y
    
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    call x%GuardTemp()
    select type(x)
    class is (par_scalar_array_t)
       select type(y)
       class is(par_scalar_array_t)
          call par_preconditioner_dd_diagonal_apply_all_unk ( op, x, y )
       end select
    end select
    call x%CleanTemp()
  end subroutine par_preconditioner_dd_diagonal_apply
  
  function par_preconditioner_dd_diagonal_is_linear(op)
    implicit none
    class(par_preconditioner_dd_diagonal_t), intent(in) :: op
    logical :: par_preconditioner_dd_diagonal_is_linear
    par_preconditioner_dd_diagonal_is_linear = .true.
  end function par_preconditioner_dd_diagonal_is_linear
  
  subroutine par_preconditioner_dd_diagonal_free(this)
    implicit none
    class(par_preconditioner_dd_diagonal_t), intent(inout) :: this
    call par_preconditioner_dd_diagonal_free_in_stages(this,free_values)
    call par_preconditioner_dd_diagonal_free_in_stages(this,free_clean)
  end subroutine par_preconditioner_dd_diagonal_free
  
end module par_preconditioner_dd_diagonal_names
