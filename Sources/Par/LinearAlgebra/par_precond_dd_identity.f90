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
module par_precond_dd_identity_names
  ! Serial modules
  use types
  use memor
  use fem_vector_names
  use fem_precond_names, only: invert_diagonal, apply_diagonal, extract_diagonal

  ! Parallel modules
  use par_vector_names
  use par_matrix_names
  use par_context_names
  use par_environment_names
  use dof_distribution_names
  use psb_penv_mod

  ! Abstract modules
  use base_operand_names
  use base_operator_names

# include "debug.i90"
  
  implicit none
  private

  type, extends(base_operator) :: par_precond_dd_identity
     ! Reference to parallel matrix
     type( par_matrix ), pointer     :: p_mat => NULL()   
     real(rp)          , allocatable :: d(:)            ! Inverse of main diagonal
   contains
     procedure :: apply     => par_precond_dd_identity_apply_tbp
     procedure :: apply_fun => par_precond_dd_identity_apply_fun_tbp
     procedure :: free      => par_precond_dd_identity_free_tbp
  end type par_precond_dd_identity
  
  ! Types
  public :: par_precond_dd_identity

  ! Functions
  public :: par_precond_dd_identity_create, par_precond_dd_identity_free, &
            par_precond_dd_identity_ass_struct, par_precond_dd_identity_fill_val, &
            par_precond_dd_identity_apply_all_unk

  contains

  !=============================================================================
  subroutine par_precond_dd_identity_create (p_matrix, p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_matrix)             , target, intent(in)  :: p_matrix
    type(par_precond_dd_identity)        , intent(out) :: p_prec_dd_identity

    
    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist) )

    p_prec_dd_identity%p_mat    => p_matrix

  end subroutine par_precond_dd_identity_create

  !=============================================================================
  subroutine par_precond_dd_identity_ass_struct (p_matrix, p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_matrix)             , target, intent(in)    :: p_matrix
    type(par_precond_dd_identity)        , intent(inout) :: p_prec_dd_identity

    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist) )

    p_prec_dd_identity%p_mat    => p_matrix

  end subroutine par_precond_dd_identity_ass_struct

  !=============================================================================
  subroutine par_precond_dd_identity_fill_val (p_matrix, p_prec_dd_identity)
    implicit none
    ! Parameters
    type(par_matrix)             , target, intent(in)    :: p_matrix
    type(par_precond_dd_identity), target, intent(inout) :: p_prec_dd_identity

    ! Locals
    integer(ip)       :: neq
    type (par_vector) :: p_vec

    assert ( associated(p_matrix%p_env) )
    assert ( p_matrix%p_env%created )
    assert ( associated(p_matrix%dof_dist) )

    p_prec_dd_identity%p_mat  => p_matrix

  end subroutine par_precond_dd_identity_fill_val

  !=============================================================================
  subroutine par_precond_dd_identity_apply_all_unk (p_prec_dd_identity, x, y)
    implicit none
    ! Parameters
    type(par_precond_dd_identity) , intent(in)    :: p_prec_dd_identity
    type(par_vector)              , intent(in)    :: x
    type(par_vector)              , intent(inout) :: y

    assert ( associated(p_prec_dd_identity%p_mat) )
    assert ( associated(p_prec_dd_identity%p_mat%p_env) )
    assert ( p_prec_dd_identity%p_mat%p_env%created )
    assert ( associated(p_prec_dd_identity%p_mat%dof_dist) )

    if(p_prec_dd_identity%p_mat%p_env%p_context%iam<0) return

    ! Comm
    if ( x%state == part_summed ) then
       call y%comm()
    end if

  end subroutine par_precond_dd_identity_apply_all_unk

  !=============================================================================
  subroutine  par_precond_dd_identity_free (p_prec_dd_identity, mode)
    implicit none
    ! Parameters
    type(par_precond_dd_identity),  intent(inout) :: p_prec_dd_identity
    integer(ip)                  ,  intent(in)    :: mode

    assert ( associated(p_prec_dd_identity%p_mat) )
    assert ( associated(p_prec_dd_identity%p_mat%p_env) )
    assert ( p_prec_dd_identity%p_mat%p_env%created )
    assert ( associated(p_prec_dd_identity%p_mat%dof_dist) )

    nullify ( p_prec_dd_identity%p_mat )

  end subroutine par_precond_dd_identity_free

  !=============================================================================
  subroutine par_precond_dd_identity_apply_tbp (op, x, y)
    implicit none
    ! Parameters
    class(par_precond_dd_identity)    , intent(in)    :: op
    class(base_operand)   , intent(in)    :: x
    class(base_operand)   , intent(inout) :: y
        
    call x%GuardTemp()
    
    select type(x)
    class is (par_vector)
       select type(y)
       class is(par_vector)
          call par_precond_dd_identity_apply_all_unk ( op, x, y )
       class default
          write(0,'(a)') 'fem_matrix%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'par_precond_dd_identity%apply: unsupported x class'
       check(1==0)
    end select
    
    call x%CleanTemp()
  end subroutine par_precond_dd_identity_apply_tbp
  
  
  !=============================================================================
  function par_precond_dd_identity_apply_fun_tbp (op, x) result(y)
    implicit none
    ! Parameters
    class(par_precond_dd_identity), intent(in)   :: op
    class(base_operand), intent(in)  :: x
    class(base_operand), allocatable :: y
    type(par_vector), allocatable :: local_y
    
    call x%GuardTemp()
    
    select type(x)
    class is (par_vector)
       allocate(local_y)
       call par_vector_alloc ( x%dof_dist, x%p_env, local_y)
       call par_precond_dd_identity_apply_all_unk ( op, x, local_y )
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'par_precond_dd_identity%apply_fun: unsupported x class'
       check(1==0)
    end select
    
    call x%CleanTemp()
  end function par_precond_dd_identity_apply_fun_tbp
  
  subroutine par_precond_dd_identity_free_tbp(this)
    implicit none
    class(par_precond_dd_identity), intent(inout) :: this
  end subroutine par_precond_dd_identity_free_tbp
  
end module par_precond_dd_identity_names
