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
module fem_vector_krylov_basis_names
  use types
  use memor
#ifdef ENABLE_BLAS
  use blas77_interfaces 
#endif
  use fem_vector_names
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

  ! fem_vector_krylov_basis
  type fem_vector_krylov_basis
     integer(ip)                :: &
        nd  = 0,                   &  ! Number of degrees of freedom, ndof
        neq = 0,                   &  ! Number of equations
        k   = 0                       ! Number of Krylov basis vectors

     integer(ip)                :: &  ! Storage layout (blk: block; scal: scalar)
        storage = blk

     real(rp), allocatable      :: &
        b(:,:,:) 
  end type fem_vector_krylov_basis

  ! Types
  public :: fem_vector_krylov_basis

  ! Functions
  public :: fem_vector_krylov_basis_alloc, fem_vector_krylov_basis_free,            & 
            fem_vector_krylov_basis_extract_view, fem_vector_krylov_basis_multidot, & 
            fem_vector_krylov_basis_multiaxpy

contains

  !=============================================================================
  subroutine fem_vector_krylov_basis_alloc (k, f_v, Q) 
    implicit none
    integer(ip)                  , intent(in)          :: k
    type(fem_vector)             , intent(in) , target :: f_v
    type(fem_vector_krylov_basis), intent(out)         :: Q 

    Q%nd      = f_v%nd
    Q%neq     = f_v%neq
    Q%k       = k
    Q%storage = f_v%storage

    if ( Q%storage == blk ) then 
      call memalloc(Q%nd, Q%neq, Q%k, Q%b, __FILE__,__LINE__)
    else if ( Q%storage == scal ) then
      call memalloc(1, Q%nd*Q%neq, Q%k, Q%b, __FILE__,__LINE__)
    end if 

    Q%b = 0.0_rp 
  end subroutine fem_vector_krylov_basis_alloc

  !=============================================================================
  subroutine fem_vector_krylov_basis_free (Q)
     implicit none
     type(fem_vector_krylov_basis), intent(inout) :: Q
    
     Q%nd      = 0 
     Q%neq     = 0
     Q%k       = 0
     Q%storage = undef_sto
     call memfree(Q%b,__FILE__,__LINE__)

  end subroutine fem_vector_krylov_basis_free

  !=============================================================================
  subroutine fem_vector_krylov_basis_extract_view (i, Q, f_v)
     implicit none
     integer(ip)     , intent(in)                      :: i
     type(fem_vector_krylov_basis), intent(in), target :: Q
     type(fem_vector), intent(out)                     :: f_v

     assert ( i >= 1 .and. i <= Q%k )

     ! fill f_v members
     f_v%nd      =  Q%nd
     f_v%neq     =  Q%neq
     f_v%mode    =  reference
     f_v%storage =  Q%storage
     f_v%b       => Q%b(:,:,i)
  end subroutine fem_vector_krylov_basis_extract_view

  !=============================================================================
  ! s <- Q_k^T * f_v, with Q_k = (Q(1), Q(2), .. Q(k))
  subroutine fem_vector_krylov_basis_multidot (k, Q, f_v, s)
     implicit none
     integer(ip)                  , intent(in) :: k
     type(fem_vector_krylov_basis), intent(in) :: Q
     type(fem_vector)             , intent(in) :: f_v
     real(rp), intent(out)                     :: s(k)

     assert ( f_v%nd  == Q%nd  )
     assert ( f_v%neq == Q%neq )
     assert ( f_v%storage == Q%storage )

#ifdef ENABLE_BLAS 
     call DGEMV(  'T', Q%neq * Q%nd, k, 1.0_rp, Q%b, &
               &  Q%neq * Q%nd, f_v%b, 1, 0.0_rp, s, 1)    
#else
     write (0,*) 'Error: fem_vector_krylov_basis_multidot was not compiled with -DENABLE_BLAS.'
     write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
     stop
#endif 
     
  end subroutine fem_vector_krylov_basis_multidot

  !=============================================================================
  ! f_v <- f_v + alpha*Q_k * s
  subroutine fem_vector_krylov_basis_multiaxpy (k, alpha, Q, s, f_v)
     implicit none
     integer(ip)                  , intent(in)    :: k
     real(rp)                     , intent(in)    :: alpha
     type(fem_vector_krylov_basis), intent(in)    :: Q
     real(rp)                     , intent(in)    :: s(k)
     type(fem_vector)             , intent(inout) :: f_v
     
     assert ( f_v%nd  == Q%nd  )
     assert ( f_v%neq == Q%neq )
     assert ( f_v%storage == Q%storage )    
 
#ifdef ENABLE_BLAS
      call DGEMV(  'N', Q%neq * Q%nd, k, alpha, Q%b, &
              &  Q%neq * Q%nd, s, 1, 1.0_rp, f_v%b, 1)
#else
     write (0,*) 'Error: fem_vector_krylov_basis_multidot was not compiled with -DENABLE_BLAS.'
     write (0,*) 'Error: You must activate this cpp macro in order to use the BLAS'
     stop    
#endif
   
  end subroutine fem_vector_krylov_basis_multiaxpy
 

end module fem_vector_krylov_basis_names
