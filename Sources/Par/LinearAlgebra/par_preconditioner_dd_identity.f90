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
  use vector_names
  use preconditioner_names, only: invert_diagonal, apply_diagonal, extract_diagonal

  ! Parallel modules
  use par_vector_names
  use par_matrix_names
  use par_context_names
  use par_environment_names
  use dof_distribution_names
use psb_penv_mod_names

  ! Abstract modules
  use base_operand_names
  use abstract_operator_names

# include "debug.i90"
  
  implicit none
  private

  type, extends(abstract_operator_t) :: par_preconditioner_dd_identity_t
     ! Reference to parallel matrix
     type( par_matrix_t ), pointer     :: p_mat => NULL()   
     real(rp)          , allocatable :: d(:)            ! Inverse of main diagonal
   contains
     procedure :: apply     => par_preconditioner_dd_identity_apply_tbp
     procedure :: apply_fun => par_preconditioner_dd_identity_apply_fun_tbp
     procedure :: free      => par_preconditioner_dd_identity_free_tbp
  end type par_preconditioner_dd_identity_t
  
  ! Types
  public :: par_preconditioner_dd_identity_t

  ! Functions
  public :: par_preconditioner_dd_identity_create, par_preconditioner_dd_identity_free, &
            par_preconditioner_dd_identity_ass_struct, par_preconditioner_dd_identity_fill_val, &
            par_preconditioner_dd_identity_apply_all_unk

  contains

  !=============================================================================
  subroutine par_preconditioner_dd_identity_create (p_matrix, p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_matrix_t)             , target, intent(in)  :: p_matrix
    type(par_preconditioner_dd_identity_t)        , intent(out) :: p_prec_dd_identity

    
    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist) )

    p_prec_dd_identity%p_mat    => p_matrix

  end subroutine par_preconditioner_dd_identity_create

  !=============================================================================
  subroutine par_preconditioner_dd_identity_ass_struct (p_matrix, p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_matrix_t)             , target, intent(in)    :: p_matrix
    type(par_preconditioner_dd_identity_t)        , intent(inout) :: p_prec_dd_identity

    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist) )

    p_prec_dd_identity%p_mat    => p_matrix

  end subroutine par_preconditioner_dd_identity_ass_struct

  !=============================================================================
  subroutine par_preconditioner_dd_identity_fill_val (p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_identity_t), target, intent(inout) :: p_prec_dd_identity

    ! Locals
    type(par_matrix_t), pointer :: p_matrix
    integer(ip)       :: neq
    type (par_vector_t) :: p_vec

    p_matrix  => p_prec_dd_identity%p_mat

    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist) )


  end subroutine par_preconditioner_dd_identity_fill_val

  !=============================================================================
  subroutine par_preconditioner_dd_identity_apply_all_unk (p_prec_dd_identity, x, y)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_identity_t) , intent(in)    :: p_prec_dd_identity
    type(par_vector_t)              , intent(in)    :: x
    type(par_vector_t)              , intent(inout) :: y

    assert ( associated(p_prec_dd_identity%p_mat) )
    assert ( associated(p_prec_dd_identity%p_mat%p_env) )
    assert ( p_prec_dd_identity%p_mat%p_env%created )
    assert ( associated(p_prec_dd_identity%p_mat%dof_dist) )

    if(p_prec_dd_identity%p_mat%p_env%p_context%iam<0) return

    call y%copy(x)
    
    ! Comm
    if ( x%state == part_summed ) then
       call y%comm()
    end if

  end subroutine par_preconditioner_dd_identity_apply_all_unk

  !=============================================================================
  subroutine  par_preconditioner_dd_identity_free (p_prec_dd_identity, mode)
    implicit none
    ! Parameters
    type(par_preconditioner_dd_identity_t),  intent(inout) :: p_prec_dd_identity
    integer(ip)                  ,  intent(in)    :: mode

    assert ( associated(p_prec_dd_identity%p_mat) )
    assert ( associated(p_prec_dd_identity%p_mat%p_env) )
    assert ( p_prec_dd_identity%p_mat%p_env%created )
    assert ( associated(p_prec_dd_identity%p_mat%dof_dist) )

    if ( mode == free_clean ) then
       nullify ( p_prec_dd_identity%p_mat )
    end if

  end subroutine par_preconditioner_dd_identity_free

  !=============================================================================
  subroutine par_preconditioner_dd_identity_apply_tbp (op, x, y)
    implicit none
    ! Parameters
    class(par_preconditioner_dd_identity_t)    , intent(in)    :: op
    class(base_operand_t)   , intent(in)    :: x
    class(base_operand_t)   , intent(inout) :: y
        
    call x%GuardTemp()
    
    select type(x)
    class is (par_vector_t)
       select type(y)
       class is(par_vector_t)
          call par_preconditioner_dd_identity_apply_all_unk ( op, x, y )
       class default
          write(0,'(a)') 'matrix_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'par_preconditioner_dd_identity_t%apply: unsupported x class'
       check(1==0)
    end select
    
    call x%CleanTemp()
  end subroutine par_preconditioner_dd_identity_apply_tbp
  
  
  !=============================================================================
  function par_preconditioner_dd_identity_apply_fun_tbp (op, x) result(y)
    implicit none
    ! Parameters
    class(par_preconditioner_dd_identity_t), intent(in)   :: op
    class(base_operand_t), intent(in)  :: x
    class(base_operand_t), allocatable :: y
    type(par_vector_t), allocatable :: local_y
    
    call x%GuardTemp()
    
    select type(x)
    class is (par_vector_t)
       allocate(local_y)
       call par_vector_alloc ( x%dof_dist, x%p_env, local_y)
       call par_preconditioner_dd_identity_apply_all_unk ( op, x, local_y )
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'par_preconditioner_dd_identity_t%apply_fun: unsupported x class'
       check(1==0)
    end select
    
    call x%CleanTemp()
  end function par_preconditioner_dd_identity_apply_fun_tbp
  
  subroutine par_preconditioner_dd_identity_free_tbp(this)
    implicit none
    class(par_preconditioner_dd_identity_t), intent(inout) :: this
  end subroutine par_preconditioner_dd_identity_free_tbp
  
end module par_preconditioner_dd_identity_names
