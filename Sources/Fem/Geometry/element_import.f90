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
module fem_element_import_class
  use types
  use memor
  implicit none
  private

  ! Element_Import
  ! Host for data needed to perform nearest neighbour communications for
  ! the distributed dual graph associated to the finite element mesh 
  type fem_element_import
     integer(ip)       ::          &
        ipart,                     &    ! Part identifier
        nparts                          ! Number of parts

     integer(ip)              :: nelem  ! Number of local elements in this part
     integer(ip)              :: nghost ! Number of ghost elements of this part

     integer(ip)              :: npadj    ! Number of adjacent parts
     integer(ip), allocatable :: lpadj(:) ! List of adjacent parts 

     integer(ip), allocatable  :: &           
        rcv_ptrs(:)                     ! How many elements does the part receive from each neighbour ?
     integer(igp), allocatable :: &           
        rcv_geids(:)                    ! Global element IDs of the elements to be received by this part  

     integer, allocatable ::     &           
        snd_ptrs(:)                     ! How many elements does does the part send to each neighbour?
     integer(ip), allocatable :: & 
        snd_leids(:)                    ! Local elements IDs of the elements to be sent to each neighbour ?
  end type fem_element_import

  ! Types
  public :: fem_element_import

  ! Functions
  public :: fem_element_import_free, fem_element_import_print

contains

  !=============================================================================
  subroutine fem_element_import_free ( element_import )
    !-----------------------------------------------------------------------
    ! This routine frees an import object
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_element_import), intent(inout)  :: element_import

    ! Free lpadj vector
    call memfree ( element_import%lpadj,__FILE__,__LINE__)

    ! Free snd/rcv info.
    call memfree ( element_import%rcv_ptrs,__FILE__,__LINE__)
    call memfree ( element_import%rcv_geids,__FILE__,__LINE__)

    call memfree ( element_import%snd_ptrs,__FILE__,__LINE__)
    call memfree ( element_import%snd_leids,__FILE__,__LINE__)
 
    return

  end subroutine fem_element_import_free

  !=============================================================================
  subroutine fem_element_import_print (lu_out, element_import)
    !-----------------------------------------------------------------------
    ! This routine prints an element_import object
    !-----------------------------------------------------------------------
    use types
    implicit none


    ! Parameters
    type(fem_element_import)   , intent(in)  :: element_import
    integer(ip)        , intent(in)  :: lu_out

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then

       write(lu_out,'(a)') '*** begin fem_element_import data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
            &  element_import%nparts

       write(lu_out,'(a,i10)') 'Number of neighbours:', &
            &  element_import%npadj

       write(lu_out,'(a)') 'List of neighbours:'
       write(lu_out,'(10i10)') element_import%lpadj(1:element_import%npadj)

       write(lu_out,'(a)') 'Rcv_ptrs:'
       write(lu_out,'(10i10)') element_import%rcv_ptrs(1:element_import%npadj+1)
       write(lu_out,'(a)') 'Rcv_geids:'
       write(lu_out,'(10i10)') element_import%rcv_geids
       

       write(lu_out,'(a)') 'Snd_ptrs:'
       write(lu_out,'(10i10)') element_import%snd_ptrs(1:element_import%npadj+1)
       write(lu_out,'(a)') 'Snd_leids:'
       write(lu_out,'(10i10)') element_import%snd_leids

    end if

    return

  end subroutine fem_element_import_print

end module fem_element_import_class
