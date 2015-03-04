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
module par_block_vector_names
  ! Serial modules
  use types
  use memor
  use fem_block_vector_names
  use fem_vector_names
    
  ! Parallel modules
  use par_vector_names

  implicit none
# include "debug.i90"

  private
 
  !=============================================================
  ! TO-CONSIDER:
  ! 
  ! x par_block_vector_dot and par_block_vector_nrm2 are not 
  !   optimal as a global communication (i.e., allreduce) is 
  !   required for each block. A more clever strategy could 
  !   replicate the body of par_vector_dot on each block of 
  !   par_block_vector_dot, and reduce the number of global 
  !   communications to one. However, with this strategy the code 
  !   would become significantly dirty. Does it really pay off ?
  ! 
  ! x Support for this parallel data structure in integrate.i90 
  !   would eliminate f_blk_vector member of par_block_vector and 
  !   par_block_vector_fill_complete method.
  !=============================================================

  !=============================================================
  ! TO-DO (CRUCIAL):
  !
  ! x This implementation of par_block_vector DOES NOT WORK for
  !   element-based data distributions, as the following 
  !   operations of par_vector:
  !     * par_vector_comm
  !     * par_vector_weight
  !     * par_vector_nrm2 (nrm2_part_summed, nrm2_full_summed)
  !     * par_vector_dot  (weighted_dot)
  !   are well-defined if and only if the I/O par_vectors
  !   are vectors on the interface. NOTE that the blocks of 
  !   par_block_vector are not vectors on the interface 
  !   (i.e., they are GLOBAL). 
  ! 
  !=============================================================


  ! par_vector
  type par_block_vector
     integer(ip)                   :: nblocks = 0
     type(par_vector), allocatable :: blocks(:)

     ! **IMPORTANT NOTE**: This is an auxiliary data 
     ! structure provided in order to use SERIAL 
     ! fem_block_vector assembly routines. The blocks of this 
     ! data structure are just VIEWS to the corresponding 
     ! counterparts in type(par_vector), allocatable :: blocks(:).
     ! This is required because currently integrate.i90 only
     ! accepts fem* data structures. If we provided support for 
     ! par* data structures in integrate.i90 we would not require 
     ! this aux. data structure
     type(fem_block_vector)        :: f_blk_vector
     logical                       :: fill_completed
  end type par_block_vector

  ! Types
  public :: par_block_vector

  ! Functions
  public :: par_block_vector_free,          par_block_vector_alloc,             &
            par_block_vector_fill_complete, par_block_vector_create_view,       &
            par_block_vector_clone,         par_block_vector_comm,              &
            par_block_vector_weight,                                            &
            par_block_vector_get_external_data,                                 &
            par_block_vector_dot,           par_block_vector_nrm2,              &
            par_block_vector_copy,          par_block_vector_zero,              &
            par_block_vector_init,          par_block_vector_scale,             & 
            par_block_vector_mxpy,          par_block_vector_axpy,              &
            par_block_vector_aypx,          par_block_vector_pxpy,              &
            par_block_vector_pxmy,          par_block_vector_print     
contains

  !=============================================================================
  subroutine par_block_vector_free (bvec)
    implicit none
    type(par_block_vector), intent(inout) :: bvec

    ! Locals
    integer(ip) :: ib

    do ib=1, bvec%nblocks
       call par_vector_free ( bvec%blocks(ib) )
    end do

    bvec%nblocks = 0
    deallocate( bvec%blocks )
    if ( bvec%fill_completed ) then
      call fem_block_vector_free ( bvec%f_blk_vector )
      bvec%fill_completed = .false.
    end if
  end subroutine par_block_vector_free

  !=============================================================================
  subroutine par_block_vector_alloc (nblocks, bvec)
    implicit none
    integer(ip)           , intent(in)  :: nblocks
    type(par_block_vector), intent(out) :: bvec

    bvec%nblocks        = nblocks
    bvec%fill_completed = .false.
    allocate ( bvec%blocks(nblocks) )
  end subroutine par_block_vector_alloc

  !=============================================================================
  subroutine par_block_vector_fill_complete (bvec)
    implicit none
    ! Parameters
    type(par_block_vector), intent(inout) :: bvec
    
    ! Locals
    integer(ip) :: ib  
  
    assert ( .not. bvec%fill_completed )
  
    call fem_block_vector_alloc ( bvec%nblocks, bvec%f_blk_vector )

    do ib=1, bvec%nblocks
       call fem_vector_create_view ( bvec%blocks(ib)%f_vector,        &  
                                   & 1, bvec%blocks(ib)%f_vector%neq, &
                                   & bvec%f_blk_vector%blocks(ib))
    end do
    bvec%fill_completed = .true.
  end subroutine par_block_vector_fill_complete

  !=============================================================================
  subroutine par_block_vector_create_view (svec, start, end, tvec)
    implicit none
    ! Parameters
    type(par_block_vector), intent(in)  :: svec
    integer(ip)     , intent(in)        :: start
    integer(ip)     , intent(in)        :: end
    type(par_block_vector), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib

    call par_block_vector_alloc ( svec%nblocks, tvec )

    do ib=1, svec%nblocks
       call par_vector_create_view (svec%blocks(ib), start, end, tvec%blocks(ib))
    end do
  end subroutine par_block_vector_create_view

  !=============================================================================
  subroutine par_block_vector_clone ( svec, tvec )
    implicit none
    ! Parameters
    type(par_block_vector), intent( in ) :: svec
    type(par_block_vector), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib

    call par_block_vector_alloc ( svec%nblocks, tvec )

    do ib=1, svec%nblocks
       call par_vector_clone (svec%blocks(ib), tvec%blocks(ib))
    end do

  end subroutine par_block_vector_clone

  subroutine par_block_vector_comm ( p_vec, alpha, mode )
    implicit none

    ! Parameters
    type(par_block_vector), intent(inout)         :: p_vec
    real(rp)              , intent(in), optional  :: alpha
    integer(ip)           , intent(in), optional  :: mode

    ! Local variables
    integer(ip) :: ib
    do ib=1, p_vec%nblocks
       call par_vector_comm ( p_vec%blocks(ib), alpha, mode )
    end do

  end subroutine par_block_vector_comm

  subroutine par_block_vector_weight ( p_vec )
    implicit none

    ! Parameters
    type(par_block_vector), intent(inout)         :: p_vec

    ! Local variables
    integer(ip) :: ib
    do ib=1, p_vec%nblocks
       call par_vector_weight ( p_vec%blocks(ib) )
    end do

  end subroutine par_block_vector_weight

  subroutine par_block_vector_get_external_data (bvec)
    implicit none
    ! Parameters
    type(par_block_vector), intent(inout)  :: bvec  

    ! Locals
    integer(ip)                            :: ib

    do ib=1,bvec%nblocks
      call par_vector_get_external_data ( bvec%blocks(ib) )
    end do 

  end subroutine par_block_vector_get_external_data

  !=============================================================================
  subroutine par_block_vector_dot (x, y, t)
    implicit none
    ! Parameters
    type(par_block_vector), intent(in)  :: x
    type(par_block_vector), intent(in)  :: y
    real(rp)              , intent(out) :: t
     
    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )

    t = 0.0_rp
    do ib=1,x%nblocks
      call par_vector_dot ( x%blocks(ib), y%blocks(ib), aux )
      t = t + aux
    end do 
  end subroutine par_block_vector_dot

  !=============================================================================
  subroutine par_block_vector_nrm2(x,t)
    implicit none
    type(par_block_vector), intent(in)  :: x
    real(rp)              , intent(out) :: t

    ! p_part%p_context is required within this subroutine
    assert ( associated(x%blocks(1)%p_part%p_context) )
    assert ( x%blocks(1)%p_part%p_context%created .eqv. .true.)

    if(x%blocks(1)%p_part%p_context%iam<0) return

    call par_block_vector_dot (x, x, t)
    t = sqrt(t)
  end subroutine par_block_vector_nrm2

  !=============================================================================
  subroutine par_block_vector_copy(x,y)
    implicit none
    type(par_block_vector), intent(in)    :: x
    type(par_block_vector), intent(inout) :: y

    ! Locals
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, x%nblocks
      call par_vector_copy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_copy

  subroutine par_block_vector_zero(y)
    implicit none
    type(par_block_vector), intent(inout) :: y
    ! Locals
    integer(ip) :: ib

    do ib=1, y%nblocks
      call par_vector_zero ( y%blocks(ib) )
    end do 

  end subroutine par_block_vector_zero

  subroutine par_block_vector_init(alpha, y)
    implicit none
    type(par_block_vector), intent(inout) :: y 
    real(rp), intent(in)                  :: alpha  
    ! Locals
    integer(ip)                           :: ib

    do ib=1, y%nblocks
      call par_vector_init ( alpha, y%blocks(ib) )
    end do    
  end subroutine par_block_vector_init
  
  subroutine par_block_vector_scale(t, x, y)
    implicit none
    ! Parameters 
    real(rp)              , intent(in)    :: t
    type(par_block_vector), intent(in)    :: x
    type(par_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_scale ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine par_block_vector_scale

  subroutine par_block_vector_mxpy(x,y)
    implicit none
    type(par_block_vector), intent(in)    :: x
    type(par_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_mxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_mxpy
  subroutine par_block_vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(par_block_vector), intent(in)    :: x
    type(par_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_axpy ( t, x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_axpy

  subroutine par_block_vector_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(par_block_vector), intent(in)    :: x
    type(par_block_vector), intent(inout) :: y

    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_aypx ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine par_block_vector_aypx

  subroutine par_block_vector_pxpy(x,y)
    implicit none
    type(par_block_vector), intent(in)    :: x
    type(par_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_pxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_pxpy

  subroutine par_block_vector_pxmy(x,y)
    implicit none
    type(par_block_vector), intent(in)    :: x
    type(par_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib
       
    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_pxmy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_pxmy

  subroutine par_block_vector_print (luout, x)
    implicit none
    type(par_block_vector), intent(in) :: x
    integer(ip)           , intent(in) :: luout
    
    ! Locals
    integer(ip) :: ib

    do ib=1, x%nblocks
       write (*,*) 'Block-vector ', ib
       call par_vector_print ( luout, x%blocks(ib) )
    end do 
  end subroutine par_block_vector_print

end module par_block_vector_names
