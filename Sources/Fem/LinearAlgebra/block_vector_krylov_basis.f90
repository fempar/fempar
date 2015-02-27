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
module fem_block_vector_krylov_basis_names
  use types
  use memor
  use fem_vector_names
  use fem_block_vector_names
  use fem_vector_krylov_basis_names
  implicit none
# include "debug.i90"

  !=============================================================
  ! TODO:
  ! 
  ! x Call to BLAS double precision or single precision 
  !   subroutines depending on the value of the rp parameter. 
  !   Currently we are always calling double precision variants 
  !   of the BLAS subroutines.
  ! 
  !=============================================================

  private

  type fem_block_vector_krylov_basis
     integer(ip)                                :: nblocks = 0
     type(fem_vector_krylov_basis), allocatable :: blocks(:)
  end type fem_block_vector_krylov_basis

  ! Types
  public :: fem_block_vector_krylov_basis

  ! Functions
  public :: fem_block_vector_krylov_basis_alloc,        fem_block_vector_krylov_basis_free,     & 
            fem_block_vector_krylov_basis_extract_view, fem_block_vector_krylov_basis_multidot, & 
            fem_block_vector_krylov_basis_multiaxpy

contains

  !=============================================================================
  subroutine fem_block_vector_krylov_basis_alloc (k, f_v, Q)
    implicit none
    ! Parameters
    integer (ip)                        , intent(in)  :: k
    type (fem_block_vector)             , intent(in)  :: f_v
    type (fem_block_vector_krylov_basis), intent(out) :: Q

    ! Locals
    integer(ip) :: ib 

    Q%nblocks = f_v%nblocks
    allocate ( Q%blocks (Q%nblocks) )
    do ib=1, f_v%nblocks
       call fem_vector_krylov_basis_alloc ( k, f_v%blocks(ib), Q%blocks(ib) )
    end do
    
  end subroutine fem_block_vector_krylov_basis_alloc

  !=============================================================================
  subroutine fem_block_vector_krylov_basis_free (Q)
     implicit none
     ! Parameters
     type(fem_block_vector_krylov_basis), intent(inout) :: Q

     ! Locals
     integer(ip) :: ib

     do ib=1, Q%nblocks
       call fem_vector_krylov_basis_free ( Q%blocks(ib) )
     end do

     deallocate ( Q%blocks )
  end subroutine fem_block_vector_krylov_basis_free

  !=============================================================================
  subroutine fem_block_vector_krylov_basis_extract_view (i, Q, f_v)
     implicit none
     integer(ip)     , intent(in)                      :: i
     type(fem_block_vector_krylov_basis), intent(in), target :: Q
     type(fem_block_vector), intent(out)                     :: f_v
     ! Locals
     integer(ip) :: ib

    call fem_block_vector_alloc ( Q%nblocks, f_v )

     do ib=1, f_v%nblocks
       call fem_vector_krylov_basis_extract_view ( i, Q%blocks(ib), f_v%blocks(ib) )
     end do
  end subroutine fem_block_vector_krylov_basis_extract_view

  !=============================================================================
  ! s <- Q_k^T * f_v, with Q_k = (Q(1), Q(2), .. Q(k))
  subroutine fem_block_vector_krylov_basis_multidot (k, Q, f_v, s)
     implicit none
     ! Parameters
     integer(ip)                        , intent(in) :: k
     type(fem_block_vector_krylov_basis), intent(in) :: Q
     type(fem_block_vector)             , intent(in) :: f_v
     real(rp), intent(out)                           :: s(k)
 
     ! Locals 
     real(rp)    :: p(k)
     integer(ip) :: ib

     s = 0.0_rp
     do ib=1, f_v%nblocks
       call fem_vector_krylov_basis_multidot ( k, Q%blocks(ib), f_v%blocks(ib), p )
       s = s + p 
     end do
  end subroutine fem_block_vector_krylov_basis_multidot

  !=============================================================================
  ! f_v <- f_v + alpha*Q_k * s
  subroutine fem_block_vector_krylov_basis_multiaxpy (k, alpha, Q, s, f_v)
     implicit none
     integer(ip)                        , intent(in)    :: k
     real(rp)                           , intent(in)    :: alpha
     type(fem_block_vector_krylov_basis), intent(in)    :: Q
     real(rp)                           , intent(in)    :: s(k)
     type(fem_block_vector)             , intent(inout) :: f_v

     ! Locals
     integer(ip) :: ib 
     do ib=1, f_v%nblocks
       call fem_vector_krylov_basis_multiaxpy ( k, alpha, Q%blocks(ib), s, f_v%blocks(ib) )
     end do
  end subroutine fem_block_vector_krylov_basis_multiaxpy
 
end module fem_block_vector_krylov_basis_names
