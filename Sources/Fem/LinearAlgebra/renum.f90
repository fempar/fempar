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
module renum_names
  use types
  use memor
  use stdio
  implicit none
  private

  type renum
     integer(ip)                :: &
        n                                    ! Size of permutation arrays
     integer(ip), allocatable ::   &
        lperm(:),                  &         ! Permutation array (lperm(old)=new)
        iperm(:)                             ! Inverse permutation array (iperm(new)=old)
  end type renum

  interface renum_apply
     module procedure renum_apply_r1, renum_apply_r2, renum_apply_i1, renum_apply_i1_igp, renum_apply_i2
  end interface

  ! Types
  public :: renum

  ! Functions
  public :: renum_alloc, renum_copy, renum_free, renum_by_sets, renum_check, renum_write, &
       &    renum_read, renum_compose_name, renum_apply, renum_identity, renum_inverse

contains
  !=============================================================================
  !=============================================================================
  subroutine renum_by_sets(nsets,lvset,ren)
    !-----------------------------------------------------------------------
    ! This routine 
    !-----------------------------------------------------------------------
    implicit none
    type(renum) , intent(inout) :: ren
    integer(ip)        , intent(in)    :: nsets
    integer(ip)        , intent(in)    :: lvset(ren%n)
    ! Local variables
    integer(ip)        , allocatable   :: iwork(:)

    call memalloc (nsets+1, iwork, __FILE__,__LINE__)
    call renum_by_sets_aux(ren%n,nsets,iwork,lvset,ren%lperm,ren%iperm)
    call memfree (iwork,__FILE__,__LINE__)

  end subroutine renum_by_sets
  !-----------------------------------------------------------------------------
  subroutine renum_by_sets_aux(np,nsets,lwork,lvset,lr,ir)
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

  end subroutine renum_by_sets_aux

  !=============================================================================
  subroutine renum_check(ren)
    implicit none
    ! Parameters
    type(renum), intent(in) :: ren

    ! Local variables
    integer(ip)               :: i
    integer(ip), allocatable  :: occurs(:)

    call memalloc ( ren%n, occurs,__FILE__,__LINE__)

    occurs   = 0
    do i = 1, ren%n
       if (ren%lperm (i) < 1 .or. ren%lperm (i) > ren%n) then
          write(*,*) '** [Fempar Warning] ** renum_check: node id ', i,      &
             &          ' is un-correctly mapped to node id', ren%lperm (i), &  
             &          ' in ren%lperm'  ! DBG: 
       end if
       occurs (ren%lperm (i)) = occurs (ren%lperm (i)) + 1
       if (occurs(ren%lperm (i)) > 1) then
          write(*,*) '** [Fempar Warning] ** renum_check: node id ', i, &
             &     ' is listed more than once in ren%lperm'  ! DBG: 
       end if
    end do

    do i = 1, ren%n
       if (occurs(i) ==0) then
          write(*,*) '** [Fempar Warning] ** renum_check: node id ', i, &
             &          ' is not listed ren%lperm'  ! DBG: 
       end if
    end do

    occurs   = 0
    do i = 1, ren%n
       if (ren%iperm (i) < 1 .or. ren%iperm (i) > ren%n) then
          write(*,*) '** [Fempar Warning] ** renum_check: node id ', i,      &
             &          ' is un-correctly mapped to node id', ren%iperm (i), &  
             &          ' in ren%iperm'  ! DBG: 
       else
          occurs (ren%iperm (i)) = occurs (ren%iperm (i)) + 1
          if (occurs(ren%iperm (i)) > 1) then
             write(*,*) '** [Fempar Warning] ** renum_check: node id ', i, &
                  &     ' is listed more than once in ren%iperm'  ! DBG: 
          end if
       end if
    end do

    do i = 1, ren%n
       if (occurs(i) ==0) then
          write(*,*) '** [Fempar Warning] ** renum_check: node id ', i, &
             &          ' is not listed ren%iperm'  ! DBG: 
       end if
    end do

    call memfree (occurs,__FILE__,__LINE__)

  end subroutine renum_check

  !=============================================================================
  subroutine renum_alloc(n,ren)
    implicit none
    integer(ip), intent(in)         :: n
    type(renum), intent(out) :: ren

    ren%n=n
    if(ren%n>0) then
       call memalloc (ren%n, ren%lperm, __FILE__,__LINE__)
       call memalloc (ren%n, ren%iperm, __FILE__,__LINE__)
       call renum_identity(ren%n,ren%lperm)
       call renum_inverse(ren%n,ren%lperm,ren%iperm)
    else
       write(*,*) 'Cannot allocate renum of',n,'components'
       stop
    end if
  end subroutine renum_alloc

  subroutine renum_copy(iren,oren)
    implicit none
    type(renum), intent(in)  :: iren
    type(renum), intent(out) :: oren

    oren%n = iren%n
    call memalloc (oren%n, oren%lperm, __FILE__,__LINE__)
    call memalloc (oren%n, oren%iperm, __FILE__,__LINE__)

    oren%lperm = iren%lperm
    oren%iperm = iren%iperm

  end subroutine renum_copy

  !-----------------------------------------------------------------------------
  subroutine renum_identity(n,l)
    implicit none
    integer(ip), intent(in)  :: n
    integer(ip), intent(out) :: l(n)
    integer(ip) :: i
    do i=1,n
       l(i)=i
    end do
  end subroutine renum_identity
  !-----------------------------------------------------------------------------
  subroutine renum_inverse(n,l1,l2)
    integer(ip), intent(in)  :: n,l1(n)
    integer(ip), intent(out) :: l2(n)
    integer(ip) :: i
    do i=1,n
       l2(l1(i))=i
    end do
  end subroutine renum_inverse
  !=============================================================================
  subroutine renum_free(ren)
    implicit none
    type(renum), intent(inout) :: ren
    if(ren%n>0) then
       call memfree (ren%lperm,__FILE__,__LINE__)
       call memfree (ren%iperm,__FILE__,__LINE__)
    end if
  end subroutine renum_free

  !=============================================================================
  subroutine renum_write(file_path,ren)
    ! Parameters
    character *(*)     , intent(in)  :: file_path
    type(renum)        , intent(in)  :: ren
    !-----------------------------------------------------------------------
    ! This routine writes a renumeration object to file file_path
    !-----------------------------------------------------------------------
    ! Locals 
    integer :: lunio
    lunio = io_open (file_path, 'write')

    write ( lunio, * ) ren%n
    write ( lunio, * ) ren%lperm
    write ( lunio, * ) ren%iperm

    call io_close (lunio)
  end subroutine renum_write

  !=============================================================================
  subroutine renum_read (file_path, ren)
    ! Parameters
    character *(*)     , intent(in)     :: file_path
    type(renum), intent(out)    :: ren
    !-----------------------------------------------------------------------
    ! This routine reads a renumeration object from file file_path
    !-----------------------------------------------------------------------
    ! Locals 
    integer :: lunio
    lunio = io_open (file_path, 'read', status='old')

    read ( lunio, * ) ren%n
    call memalloc (ren%n, ren%lperm, __FILE__,__LINE__)
    call memalloc (ren%n, ren%iperm, __FILE__,__LINE__)
    read ( lunio, * ) ren%lperm
    read ( lunio, * ) ren%iperm 

    call io_close(lunio)
  end subroutine renum_read

  !=============================================================================
  subroutine renum_compose_name ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.ren'
  end subroutine renum_compose_name


  !================================================================================================
  subroutine renum_apply_r1(ren, xlin, xlout)
    implicit none
    type(renum), intent(in)  :: ren
    real(rp)   , intent(in)  :: xlin  (ren%n)
    real(rp)   , intent(out) :: xlout (ren%n)
    integer(ip) :: i

    do i=1,ren%n
       xlout(i)=xlin(ren%iperm(i))
    end do
  end subroutine renum_apply_r1

  !================================================================================================
  subroutine renum_apply_r2(ld, ren, xlin, xlout)
    implicit none
    integer(ip), intent(in)  :: ld
    type(renum), intent(in)  :: ren
    real(rp)   , intent(in)  :: xlin  (ld,ren%n)
    real(rp)   , intent(out) :: xlout (ld,ren%n)
    integer(ip) :: i

    do i=1,ren%n
       xlout(:,i)=xlin(:,ren%iperm(i))
    end do
  end subroutine renum_apply_r2

  !================================================================================================
  subroutine renum_apply_i1(ren, xlin, xlout)
    implicit none
    type(renum), intent(in)  :: ren
    integer(ip), intent(in)  :: xlin  (ren%n)
    integer(ip), intent(out) :: xlout (ren%n)
    integer(ip)              :: i

    do i=1,ren%n
       xlout(i)=xlin(ren%iperm(i))
    end do
  end subroutine renum_apply_i1

  !================================================================================================
  subroutine renum_apply_i1_igp(ren, xlin, xlout)
     implicit none
     type(renum), intent(in)   :: ren
     integer(igp), intent(in)  :: xlin  (ren%n)
     integer(igp), intent(out) :: xlout (ren%n)
     integer(ip)               :: i

     do i=1,ren%n
        xlout(i)=xlin(ren%iperm(i))
     end do
  end subroutine renum_apply_i1_igp

  !================================================================================================
  subroutine renum_apply_i2(ld, ren, xlin, xlout)
    implicit none
    integer(ip), intent(in)  :: ld
    type(renum), intent(in)  :: ren
    integer(ip)   , intent(in)  :: xlin  (ld,ren%n)
    integer(ip)   , intent(out) :: xlout (ld,ren%n)
    integer(ip) :: i

    do i=1,ren%n
       xlout(:,i)=xlin(:,ren%iperm(i))
    end do
  end subroutine renum_apply_i2
    
end module renum_names
