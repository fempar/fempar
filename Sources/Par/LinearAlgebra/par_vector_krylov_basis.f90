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
module par_vector_krylov_basis_names
  ! Serial modules
use types_names
use memor_names
  use vector_krylov_basis_names

  ! Parallel modules
  use dof_distribution_names
  use psb_penv_mod_names
  use par_environment_names
  use par_vector_names

# include "debug.i90"
  
  implicit none
  private


  ! Distributed Krylov Basis (compatible with par_vector)
  type par_vector_krylov_basis_t
     ! Local view of ONLY those components of b
     ! corresponding to vertices owned by the processor  
     !type( epetra_multivector_t ) :: epmv

     type(vector_krylov_basis_t) :: f_basis
     
     ! Partially or fully summed
     integer(ip)  :: state = undefined

     ! DoF distribution control info.
     type ( dof_distribution_t ), pointer  :: dof_dist => NULL()
     type ( par_environment_t ) , pointer  :: p_env => NULL()
  end type par_vector_krylov_basis_t

  ! Types
  public :: par_vector_krylov_basis_t


  ! Functions
  public :: par_vector_krylov_basis_alloc, par_vector_krylov_basis_free,            & 
            par_vector_krylov_basis_extract_view, par_vector_krylov_basis_multidot, & 
            par_vector_krylov_basis_multiaxpy
contains

  !=============================================================================
  subroutine par_vector_krylov_basis_alloc (k, p_v, Q)
    implicit none
    integer(ip)     , intent(in)               :: k
    type(par_vector_t), intent(in) , target      :: p_v
    type(par_vector_krylov_basis_t), intent(out) :: Q 

    ! dof_dist and p_env%p_context is required within this subroutine
    assert ( associated(p_v%dof_dist)          )
    assert ( associated(p_v%p_env%p_context) )
    assert ( p_v%p_env%p_context%created .eqv. .true.)
    Q%dof_dist => p_v%dof_dist
    Q%p_env    => p_v%p_env


    if(p_v%p_env%p_context%iam<0) return

    call vector_krylov_basis_alloc ( k, p_v%f_vector, Q%f_basis )

    Q%state = p_v%state
  end subroutine par_vector_krylov_basis_alloc

  !=============================================================================
  subroutine par_vector_krylov_basis_free (Q)
    implicit none
    type(par_vector_krylov_basis_t), intent(inout) :: Q

    ! The routine requires the partition/context info
    assert ( associated( Q%dof_dist ) )
    assert ( associated( Q%p_env%p_context ) )
    assert ( Q%p_env%p_context%created .eqv. .true.)
    if(Q%p_env%p_context%iam<0) return

    Q%state = undefined

    call vector_krylov_basis_free ( Q%f_basis )

  end subroutine par_vector_krylov_basis_free

  !=============================================================================
  subroutine par_vector_krylov_basis_extract_view (i, Q, p_v)
    implicit none
    integer(ip)     , intent(in)                      :: i
    type(par_vector_krylov_basis_t), intent(in), target :: Q
    type(par_vector_t), intent(out)                     :: p_v

    ! The routine requires the partition/context info
    assert ( associated( Q%dof_dist ) )
    assert ( associated( Q%p_env%p_context ) )
    assert ( Q%p_env%p_context%created .eqv. .true.)

    ! 2. fill p_v members
    p_v%dof_dist => Q%dof_dist

    if(Q%p_env%p_context%iam<0) return

    ! 1. fill p_v%f_vector members
    call vector_krylov_basis_extract_view (i, Q%f_basis, p_v%f_vector )

    p_v%state  =  Q%state
  end subroutine par_vector_krylov_basis_extract_view

  !=============================================================================
  ! s <- Q_k^T * p_v, with Q_k = (Q(1), Q(2), .. Q(k))
  subroutine par_vector_krylov_basis_multidot (k, Q, p_v, s)
     implicit none
     ! Parameters 
     integer(ip), intent(in)                   :: k
     type(par_vector_krylov_basis_t), intent(in) :: Q
     type(par_vector_t)             , intent(in) :: p_v
     real(rp), intent(out)                     :: s(k)

     ! Locals
     integer(ip)                               :: i
     type(par_vector_t)                          :: p_v_w
     
     ! The routine requires the partition/context info
     assert ( associated( Q%dof_dist ) )
     assert ( associated( Q%p_env%p_context ) )
     assert ( Q%p_env%p_context%created .eqv. .true.)
     if(Q%p_env%p_context%iam<0) return
     
     assert ( p_v%state == Q%state )
     
     if ( p_v%state == full_summed ) then
        ! ******** p_v%f_vector should be weighted in advance !!!!
        call par_vector_clone  ( p_v, p_v_w )
        call par_vector_copy   ( p_v, p_v_w )
        call par_vector_weight ( p_v_w )
     else if ( p_v%state == part_summed ) then
        ! ******** p_v%f_vector should be comm. in advance !!!!
        call par_vector_clone  ( p_v, p_v_w )
        call par_vector_copy   ( p_v, p_v_w )
        call par_vector_comm   ( p_v_w )
     end if
     
     call vector_krylov_basis_multidot ( k, Q%f_basis, p_v_w%f_vector, s )
     
     call par_vector_free  ( p_v_w )
     
     ! Reduce-sum local dot products on all processes
     call psb_sum ( Q%p_env%p_context%icontxt, s )
   end subroutine par_vector_krylov_basis_multidot

  !=============================================================================
  ! p_v <- p_v + alpha * Q_k * s
  subroutine par_vector_krylov_basis_multiaxpy (k, alpha, Q, s, p_v)
     implicit none
     ! Parameters
     integer(ip), intent(in)                      :: k
     real(rp)   , intent(in)                      :: alpha
     type(par_vector_krylov_basis_t), intent(in)    :: Q
     real(rp), intent(in)                         :: s(k)
     type(par_vector_t)             , intent(inout) :: p_v

     ! Locals
     integer(ip)                                  :: i


     ! The routine requires the partition/context info
     assert ( associated( Q%dof_dist ) )
     assert ( associated( Q%p_env%p_context ) )
     assert ( Q%p_env%p_context%created .eqv. .true.)
     if(Q%p_env%p_context%iam<0) return

     assert ( p_v%state == Q%state )
     call vector_krylov_basis_multiaxpy (k, alpha, Q%f_basis, s, p_v%f_vector)

  end subroutine par_vector_krylov_basis_multiaxpy
 
end module par_vector_krylov_basis_names
