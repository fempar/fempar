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
module renumbering_names
  use types_names
  use memor_names
  use stdio_names
  implicit none
# include "debug.i90"
  private

  type renumbering_t
     integer(ip)                :: &
        n                                    ! Size of permutation arrays
     integer(ip), allocatable ::   &
        lperm(:),                  &         ! Permutation array (lperm(old)=new)
        iperm(:)                             ! Inverse permutation array (iperm(new)=old)
  end type renumbering_t

  interface renumbering_apply
     module procedure renumbering_apply_r1, renumbering_apply_r2, renumbering_apply_i1, renumbering_apply_i1_igp, renumbering_apply_i2
  end interface

  ! Types
  public :: renumbering_t

  ! Functions
  public :: renumbering_alloc, renumbering_copy, renumbering_free, renumbering_by_sets, renumbering_check, renumbering_write, &
       &    renumbering_read, renumbering_compose_name, renumbering_apply, renumbering_identity, renumbering_inverse

contains
  !=============================================================================
  !=============================================================================
  subroutine renumbering_by_sets(nsets,lvset,renumbering)
    !-----------------------------------------------------------------------
    ! This routine 
    !-----------------------------------------------------------------------
    implicit none
    type(renumbering_t) , intent(inout) :: renumbering
    integer(ip)        , intent(in)    :: nsets
    integer(ip)        , intent(in)    :: lvset(renumbering%n)
    ! Local variables
    integer(ip)        , allocatable   :: iwork(:)

    call memalloc (nsets+1, iwork, __FILE__,__LINE__)
    call renumbering_by_sets_aux(renumbering%n,nsets,iwork,lvset,renumbering%lperm,renumbering%iperm)
    call memfree (iwork,__FILE__,__LINE__)

  end subroutine renumbering_by_sets
  !-----------------------------------------------------------------------------
  subroutine renumbering_by_sets_aux(np,nsets,lwork,lvset,lr,ir)
    !-----------------------------------------------------------------------
    ! This routine 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: np, nsets
    integer(ip), intent(in)    :: lvset(np)
    integer(ip), intent(inout) :: lwork(nsets+1)
    integer(ip), intent(out)   :: lr(np), ir(np)
    ! Local variables
    integer(ip) :: i

    ! Compute how many vertices have been assigned to each set
    lwork=0
    do i=1,np
       lwork(lvset(i)+1)=lwork(lvset(i)+1)+1
    end do

    ! Compute the starting node identifier of each set 
    lwork(1) = 1
    do i=1,nsets
       lwork(i+1)=lwork(i+1)+lwork(i)
    end do

    ! Map old node identifiers to new node identifiers. lwork(lvset(i)) is 
    ! the next node identifier to be mapped within the ith-set
    do i=1,np
       lr(i)=lwork(lvset(i))
       lwork(lvset(i))=lwork(lvset(i))+1
    end do

    ! Compute inverse renumbering, i.e., map new node
    ! identifiers to old node identifiers
    do i=1,np
       ir(lr(i))=i
    end do

  end subroutine renumbering_by_sets_aux

  !=============================================================================
  subroutine renumbering_check(renumbering)
    implicit none
    ! Parameters
    type(renumbering_t), intent(in) :: renumbering

    ! Local variables
    integer(ip)               :: i
    integer(ip), allocatable  :: occurs(:)

    call memalloc ( renumbering%n, occurs,__FILE__,__LINE__)

    occurs   = 0
    do i = 1, renumbering%n
       if (renumbering%lperm (i) < 1 .or. renumbering%lperm (i) > renumbering%n) then
          write(*,*) '** [Fempar Warning] ** renumbering_check: node id ', i,      &
             &          ' is un-correctly mapped to node id', renumbering%lperm (i), &  
             &          ' in ren%lperm'  ! DBG: 
       end if
       occurs (renumbering%lperm (i)) = occurs (renumbering%lperm (i)) + 1
       if (occurs(renumbering%lperm (i)) > 1) then
          write(*,*) '** [Fempar Warning] ** renumbering_check: node id ', i, &
             &     ' is listed more than once in renumbering%lperm'  ! DBG: 
       end if
    end do

    do i = 1, renumbering%n
       if (occurs(i) ==0) then
          write(*,*) '** [Fempar Warning] ** renumbering_check: node id ', i, &
             &          ' is not listed renumbering%lperm'  ! DBG: 
       end if
    end do

    occurs   = 0
    do i = 1, renumbering%n
       if (renumbering%iperm (i) < 1 .or. renumbering%iperm (i) > renumbering%n) then
          write(*,*) '** [Fempar Warning] ** renumbering_check: node id ', i,      &
             &          ' is un-correctly mapped to node id', renumbering%iperm (i), &  
             &          ' in renumbering%iperm'  ! DBG: 
       else
          occurs (renumbering%iperm (i)) = occurs (renumbering%iperm (i)) + 1
          if (occurs(renumbering%iperm (i)) > 1) then
             write(*,*) '** [Fempar Warning] ** renumbering_check: node id ', i, &
                  &     ' is listed more than once in renumbering%iperm'  ! DBG: 
          end if
       end if
    end do

    do i = 1, renumbering%n
       if (occurs(i) ==0) then
          write(*,*) '** [Fempar Warning] ** renumbering_check: node id ', i, &
             &          ' is not listed renumbering%iperm'  ! DBG: 
       end if
    end do

    call memfree (occurs,__FILE__,__LINE__)

  end subroutine renumbering_check

  !=============================================================================
  subroutine renumbering_alloc(n,renumbering)
    implicit none
    integer(ip), intent(in)         :: n
    type(renumbering_t), intent(out) :: renumbering

    renumbering%n=n
    if(renumbering%n>0) then
       call memalloc (renumbering%n, renumbering%lperm, __FILE__,__LINE__)
       call memalloc (renumbering%n, renumbering%iperm, __FILE__,__LINE__)
       call renumbering_identity(renumbering%n,renumbering%lperm)
       call renumbering_inverse(renumbering%n,renumbering%lperm,renumbering%iperm)
    else
       write(*,*) 'Cannot allocate renumbering of',n,'components'
       check(.false.)
    end if
  end subroutine renumbering_alloc

  subroutine renumbering_copy(irenumbering,orenumbering)
    implicit none
    type(renumbering_t), intent(in)  :: irenumbering
    type(renumbering_t), intent(out) :: orenumbering

    orenumbering%n = irenumbering%n
    call memalloc (orenumbering%n, orenumbering%lperm, __FILE__,__LINE__)
    call memalloc (orenumbering%n, orenumbering%iperm, __FILE__,__LINE__)

    orenumbering%lperm = irenumbering%lperm
    orenumbering%iperm = irenumbering%iperm

  end subroutine renumbering_copy

  !-----------------------------------------------------------------------------
  subroutine renumbering_identity(n,l)
    implicit none
    integer(ip), intent(in)  :: n
    integer(ip), intent(out) :: l(n)
    integer(ip) :: i
    do i=1,n
       l(i)=i
    end do
  end subroutine renumbering_identity
  !-----------------------------------------------------------------------------
  subroutine renumbering_inverse(n,l1,l2)
    integer(ip), intent(in)  :: n,l1(n)
    integer(ip), intent(out) :: l2(n)
    integer(ip) :: i
    do i=1,n
       l2(l1(i))=i
    end do
  end subroutine renumbering_inverse
  !=============================================================================
  subroutine renumbering_free(renumbering)
    implicit none
    type(renumbering_t), intent(inout) :: renumbering
    if(renumbering%n>0) then
       call memfree (renumbering%lperm,__FILE__,__LINE__)
       call memfree (renumbering%iperm,__FILE__,__LINE__)
    end if
  end subroutine renumbering_free

  !=============================================================================
  subroutine renumbering_write(file_path,renumbering)
    ! Parameters
    character *(*)     , intent(in)  :: file_path
    type(renumbering_t)        , intent(in)  :: renumbering
    !-----------------------------------------------------------------------
    ! This routine writes a renumbering object to file file_path
    !-----------------------------------------------------------------------
    ! Locals 
    integer :: lunio
    lunio = io_open (file_path, 'write')

    write ( lunio, * ) renumbering%n
    write ( lunio, * ) renumbering%lperm
    write ( lunio, * ) renumbering%iperm

    call io_close (lunio)
  end subroutine renumbering_write

  !=============================================================================
  subroutine renumbering_read (file_path, renumbering)
    ! Parameters
    character *(*)     , intent(in)     :: file_path
    type(renumbering_t), intent(out)    :: renumbering
    !-----------------------------------------------------------------------
    ! This routine reads a renumbering object from file file_path
    !-----------------------------------------------------------------------
    ! Locals 
    integer :: lunio
    lunio = io_open (file_path, 'read', status='old')

    read ( lunio, * ) renumbering%n
    call memalloc (renumbering%n, renumbering%lperm, __FILE__,__LINE__)
    call memalloc (renumbering%n, renumbering%iperm, __FILE__,__LINE__)
    read ( lunio, * ) renumbering%lperm
    read ( lunio, * ) renumbering%iperm 

    call io_close(lunio)
  end subroutine renumbering_read

  !=============================================================================
  subroutine renumbering_compose_name ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.ren'
  end subroutine renumbering_compose_name


  !================================================================================================
  subroutine renumbering_apply_r1(renumbering, xlin, xlout)
    implicit none
    type(renumbering_t), intent(in)  :: renumbering
    real(rp)   , intent(in)  :: xlin  (renumbering%n)
    real(rp)   , intent(out) :: xlout (renumbering%n)
    integer(ip) :: i

    do i=1,renumbering%n
       xlout(i)=xlin(renumbering%iperm(i))
    end do
  end subroutine renumbering_apply_r1

  !================================================================================================
  subroutine renumbering_apply_r2(ld, renumbering, xlin, xlout)
    implicit none
    integer(ip), intent(in)  :: ld
    type(renumbering_t), intent(in)  :: renumbering
    real(rp)   , intent(in)  :: xlin  (ld,renumbering%n)
    real(rp)   , intent(out) :: xlout (ld,renumbering%n)
    integer(ip) :: i

    do i=1,renumbering%n
       xlout(:,i)=xlin(:,renumbering%iperm(i))
    end do
  end subroutine renumbering_apply_r2

  !================================================================================================
  subroutine renumbering_apply_i1(renumbering, xlin, xlout)
    implicit none
    type(renumbering_t), intent(in)  :: renumbering
    integer(ip), intent(in)  :: xlin  (renumbering%n)
    integer(ip), intent(out) :: xlout (renumbering%n)
    integer(ip)              :: i

    do i=1,renumbering%n
       xlout(i)=xlin(renumbering%iperm(i))
    end do
  end subroutine renumbering_apply_i1

  !================================================================================================
  subroutine renumbering_apply_i1_igp(renumbering, xlin, xlout)
     implicit none
     type(renumbering_t), intent(in)   :: renumbering
     integer(igp), intent(in)  :: xlin  (renumbering%n)
     integer(igp), intent(out) :: xlout (renumbering%n)
     integer(ip)               :: i

     do i=1,renumbering%n
        xlout(i)=xlin(renumbering%iperm(i))
     end do
  end subroutine renumbering_apply_i1_igp

  !================================================================================================
  subroutine renumbering_apply_i2(ld, renumbering, xlin, xlout)
    implicit none
    integer(ip), intent(in)  :: ld
    type(renumbering_t), intent(in)  :: renumbering
    integer(ip)   , intent(in)  :: xlin  (ld,renumbering%n)
    integer(ip)   , intent(out) :: xlout (ld,renumbering%n)
    integer(ip) :: i

    do i=1,renumbering%n
       xlout(:,i)=xlin(:,renumbering%iperm(i))
    end do
  end subroutine renumbering_apply_i2
    
end module renumbering_names
