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
module fem_import_names
  use types
  use memor
  implicit none
  private

  ! Import
  ! Host for data needed to perform communications for
  ! element-based data distributions (similar to epetra_import)
  type fem_import
     integer(ip)       ::          &
        ipart,                     &    ! Part identifier
        nparts                          ! Number of parts

     integer(ip), allocatable ::   &
        owner(:)                        ! Specifies for each node of the local interface
                                        ! which part is the owner
      
     integer :: num_rcv                 ! From how many neighbours does the part receive data ?
     integer, allocatable ::     &           
        list_rcv(:),             &      ! From which neighbours does the part receive data ?
        rcv_ptrs(:)                     ! How much data does the part receive from each neighbour ?
     integer(ip), allocatable :: &           
        unpack_idx(:)                   ! Where the data received from each neighbour is copied/added 
                                        ! on the local vectors of the part ?

     integer :: num_snd                 ! To how many neighbours does the part send data ? 
     integer, allocatable ::     &           
        list_snd(:),             &      ! To which neighbours does the part send data ?
        snd_ptrs(:)                     ! How much data does the part send to each neighbour?
     integer(ip), allocatable :: & 
        pack_idx(:)                     ! Where is located the data to be sent to 
                                        ! each neighbour on the local vectors of the part ?
  end type fem_import

  ! Types
  public :: fem_import

  ! Functions
  public :: fem_import_free, fem_import_print

contains

  !=============================================================================
  subroutine fem_import_free ( import )
    !-----------------------------------------------------------------------
    ! This routine frees an import object
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_import), intent(inout)  :: import

    ! Free owner vector
    call memfree ( import%owner,__FILE__,__LINE__)

    ! Free snd/rcv info.
    call memfree ( import%list_rcv,__FILE__,__LINE__)
    call memfree ( import%rcv_ptrs,__FILE__,__LINE__)
    call memfree ( import%unpack_idx,__FILE__,__LINE__)

    call memfree ( import%list_snd,__FILE__,__LINE__)
    call memfree ( import%snd_ptrs,__FILE__,__LINE__)
    call memfree ( import%pack_idx,__FILE__,__LINE__)
 
    return

  end subroutine fem_import_free

  !=============================================================================
  subroutine fem_import_print (lu_out, import)
    !-----------------------------------------------------------------------
    ! This routine prints an import object
    !-----------------------------------------------------------------------
    use types
    implicit none


    ! Parameters
    type(fem_import)   , intent(in)  :: import
    integer(ip)        , intent(in)  :: lu_out

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then

       write(lu_out,'(a)') '*** begin fem_import data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
            &  import%nparts

!!$       write(lu_out,'(a)') 'Owner of each interface point:'
!!$       write(lu_out,'(10i10)') import%owner(:)

       write(lu_out,'(a,i10)') 'Number of parts I have to receive from:', &
            &  import%num_rcv

       write(lu_out,'(a)') 'List of parts I have to receive from:'
       write(lu_out,'(10i10)') import%list_rcv(1:import%num_rcv)

       write(lu_out,'(a)') 'Rcv_ptrs:'
       write(lu_out,'(10i10)') import%rcv_ptrs(1:import%num_rcv+1)

!!$       write(lu_out,'(a)') 'Unpack:'
!!$       write(lu_out,'(10i10)') import%unpack_idx(1:import%rcv_ptrs(import%num_rcv+1)-1)
!!$
       write(lu_out,'(a,i10)') 'Number of parts I have send to:', &
            &  import%num_snd

       write(lu_out,'(a)') 'List of parts I have to send to:'
       write(lu_out,'(10i10)') import%list_snd(1:import%num_snd)

       write(lu_out,'(a)') 'Snd_ptrs::'
       write(lu_out,'(10i10)') import%snd_ptrs(1:import%num_snd+1)

!!$       write(lu_out,'(a)') 'Pack:'
!!$       write(lu_out,'(10i10)') import%pack_idx(1:import%snd_ptrs(import%num_snd+1)-1)
!!$       write(lu_out,'(a)') '*** end fem_import data structure ***'

    end if

    return

  end subroutine fem_import_print


end module fem_import_names
