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
module par_preconditioner_dd_identity_names
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

  ! Abstract modules
  use matrix_names
  use vector_names
  use operator_names
  use fe_affine_operator_names 

  implicit none
# include "debug.i90"

  private

  type, extends(operator_t) :: par_preconditioner_dd_identity_t
     ! Reference to parallel matrix
     type( par_scalar_matrix_t ), pointer     :: p_mat => NULL()   
     real(rp)          , allocatable :: d(:)            ! Inverse of main diagonal
   contains
     procedure :: is_linear => par_preconditioner_dd_identity_is_linear
     procedure :: apply     => par_preconditioner_dd_identity_apply
     procedure :: free      => par_preconditioner_dd_identity_free
  end type par_preconditioner_dd_identity_t
  
  interface par_preconditioner_dd_identity_create
    module procedure par_preconditioner_dd_identity_create_w_par_matrix, par_preconditioner_dd_identity_create_w_fe_affine_operator
  end interface
  
  ! Types
  public :: par_preconditioner_dd_identity_t

  ! Functions
  public :: par_preconditioner_dd_identity_create, par_preconditioner_dd_identity_free, &
            par_preconditioner_dd_identity_symbolic_setup, par_preconditioner_dd_identity_numerical_setup, &
            par_preconditioner_dd_identity_apply_all_unk

  contains

  !=============================================================================
  subroutine par_preconditioner_dd_identity_create_w_par_matrix (p_matrix, p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_scalar_matrix_t)             , target, intent(in)  :: p_matrix
    type(par_preconditioner_dd_identity_t)        , intent(out) :: p_prec_dd_identity
    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist_domain) )
    p_prec_dd_identity%p_mat    => p_matrix
  end subroutine par_preconditioner_dd_identity_create_w_par_matrix
  
  !=============================================================================
  subroutine par_preconditioner_dd_identity_create_w_fe_affine_operator (fe_affine_operator, p_prec_dd_identity)
    implicit none
    ! Parameters
    type(fe_affine_operator_t)            , intent(in)  :: fe_affine_operator
    type(par_preconditioner_dd_identity_t), intent(out) :: p_prec_dd_identity
    
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

    call par_preconditioner_dd_identity_create_w_par_matrix ( p_mat, p_prec_dd_identity )
  end subroutine par_preconditioner_dd_identity_create_w_fe_affine_operator

  !=============================================================================
  subroutine par_preconditioner_dd_identity_symbolic_setup (p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_identity_t), intent(inout) :: p_prec_dd_identity
  end subroutine par_preconditioner_dd_identity_symbolic_setup

  !=============================================================================
  subroutine par_preconditioner_dd_identity_numerical_setup (p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_identity_t), intent(inout) :: p_prec_dd_identity
  end subroutine par_preconditioner_dd_identity_numerical_setup

  !=============================================================================
  subroutine par_preconditioner_dd_identity_apply_all_unk (p_prec_dd_identity, x, y)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_identity_t) , intent(in)    :: p_prec_dd_identity
    type(par_scalar_array_t)              , intent(in)    :: x
    type(par_scalar_array_t)              , intent(inout) :: y

    assert ( associated(p_prec_dd_identity%p_mat) )
    assert ( associated(p_prec_dd_identity%p_mat%p_env) )
    assert ( p_prec_dd_identity%p_mat%p_env%created )
    assert ( associated(p_prec_dd_identity%p_mat%dof_dist_domain) )
    call y%copy(x)
  end subroutine par_preconditioner_dd_identity_apply_all_unk

  !=============================================================================
  subroutine  par_preconditioner_dd_identity_free_in_stages (p_prec_dd_identity, mode)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_identity_t),  intent(inout) :: p_prec_dd_identity
    integer(ip)                  ,  intent(in)    :: mode

    assert ( associated(p_prec_dd_identity%p_mat) )
    assert ( associated(p_prec_dd_identity%p_mat%p_env) )
    assert ( p_prec_dd_identity%p_mat%p_env%created )
    assert ( associated(p_prec_dd_identity%p_mat%dof_dist_domain) )

    if ( mode == free_clean ) then
       nullify ( p_prec_dd_identity%p_mat )
    end if

  end subroutine par_preconditioner_dd_identity_free_in_stages

  !=============================================================================
  subroutine par_preconditioner_dd_identity_apply (op, x, y)
    implicit none
    ! Parameters
    class(par_preconditioner_dd_identity_t)    , intent(in)    :: op
    class(vector_t)   , intent(in)    :: x
    class(vector_t)   , intent(inout) :: y
        
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    call x%GuardTemp()
    select type(x)
    class is (par_scalar_array_t)
       select type(y)
       class is(par_scalar_array_t)
          call par_preconditioner_dd_identity_apply_all_unk ( op, x, y )
       end select
    end select
    call x%CleanTemp()
  end subroutine par_preconditioner_dd_identity_apply
  
  subroutine par_preconditioner_dd_identity_free(this)
    implicit none
    class(par_preconditioner_dd_identity_t), intent(inout) :: this
    call par_preconditioner_dd_identity_free_in_stages(this,free_clean)
  end subroutine par_preconditioner_dd_identity_free
  
  function par_preconditioner_dd_identity_is_linear(op)
    implicit none
    class(par_preconditioner_dd_identity_t), intent(in) :: op
    logical :: par_preconditioner_dd_identity_is_linear
    par_preconditioner_dd_identity_is_linear = .true.
  end function par_preconditioner_dd_identity_is_linear
  
end module par_preconditioner_dd_identity_names
