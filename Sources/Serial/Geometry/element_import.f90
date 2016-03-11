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
module element_import_names
  use types_names
  use memor_names
# include "debug.i90"  
  implicit none
  private

  ! Element_Import
  ! Host for data needed to perform nearest neighbour communications for
  ! the distributed dual graph associated to the finite element mesh 
  type element_import_t
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
  contains
     procedure, non_overridable :: get_number_neighbours   => elem_import_get_number_neighbours
     procedure, non_overridable :: get_global_neighbour_id => elem_import_get_global_neighbour_id
     procedure, non_overridable :: get_local_neighbour_id  => elem_import_get_local_neighbour_id
  end type element_import_t

  ! Types
  public :: element_import_t

  ! Functions
  public :: element_import_free, element_import_print

contains

  function elem_import_get_number_neighbours ( this )
    implicit none
    class(element_import_t), intent(in) :: this
    integer(ip)                      :: elem_import_get_number_neighbours
    elem_import_get_number_neighbours = this%npadj
  end function elem_import_get_number_neighbours
  
  function elem_import_get_global_neighbour_id ( this, local_neighbour_id )
    implicit none
    class(element_import_t), intent(in) :: this
    integer(ip)         , intent(in) :: local_neighbour_id
    integer(ip)                      :: elem_import_get_global_neighbour_id
    assert ( local_neighbour_id >=1 .and. local_neighbour_id <= this%npadj )
    elem_import_get_global_neighbour_id = this%lpadj(local_neighbour_id)
  end function elem_import_get_global_neighbour_id
  
  function elem_import_get_local_neighbour_id ( this, global_neighbour_id )
    implicit none
    class(element_import_t), intent(in) :: this
    integer(ip)            , intent(in) :: global_neighbour_id
    integer(ip)                      :: elem_import_get_local_neighbour_id
        
    do elem_import_get_local_neighbour_id = 1, this%npadj
      if ( this%lpadj(elem_import_get_local_neighbour_id) == global_neighbour_id ) return
    end do
    assert ( .not. (elem_import_get_local_neighbour_id == this%npadj+1) )
  end function elem_import_get_local_neighbour_id
  
  


  !=============================================================================
  subroutine element_import_free ( element_import )
    !-----------------------------------------------------------------------
    ! This routine frees an import object
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(element_import_t), intent(inout)  :: element_import

    ! Free lpadj vector
    call memfree ( element_import%lpadj,__FILE__,__LINE__)

    ! Free snd/rcv info.
    call memfree ( element_import%rcv_ptrs,__FILE__,__LINE__)
    call memfree ( element_import%rcv_geids,__FILE__,__LINE__)

    call memfree ( element_import%snd_ptrs,__FILE__,__LINE__)
    call memfree ( element_import%snd_leids,__FILE__,__LINE__)
 
    return

  end subroutine element_import_free

  !=============================================================================
  subroutine element_import_print (lu_out, element_import)
    !-----------------------------------------------------------------------
    ! This routine prints an element_import object
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    type(element_import_t)   , intent(in)  :: element_import
    integer(ip)        , intent(in)  :: lu_out

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then
       write(lu_out,'(a)') '*** begin element_import_t data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
            &  element_import%nparts

       write(lu_out,'(a,i10)') 'Number of neighbours:', &
            &  element_import%npadj
       
       write(lu_out,'(a,i10)') 'Number of local elements:', &
            &  element_import%nelem

       write(lu_out,'(a,i10)') 'Number of ghost elements:', &
            &  element_import%nghost

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

  end subroutine element_import_print

end module element_import_names
