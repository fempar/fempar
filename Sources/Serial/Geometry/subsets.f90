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
module elem_to_subset_id_names
  use types_names
  use memor_names
  use stdio_names
  implicit none
# include "debug.i90"
  
  private

  type elem_to_subset_id_t
     private
     integer(ip)                :: &
          nelem=0                          
     integer(ip), allocatable   :: &
          elem_to_subset_id(:)
   contains
     procedure, non_overridable :: create        => elem_to_subset_id_create
     procedure, non_overridable :: free          => elem_to_subset_id_free
     procedure, non_overridable :: get_subset_id => elem_to_subset_id_get_subset_id 
  end type elem_to_subset_id_t

  ! Types
  public :: elem_to_subset_id_t

contains

  !===============================================================================================
  subroutine elem_to_subset_id_create(this,nelem)
    implicit none
    class(elem_to_subset_id_t) , intent(inout) :: this
    integer(ip)      , intent(in)    :: nelem
    call this%free()
    this%nelem=nelem
    call memalloc (this%nelem,this%elem_to_subset_id, __FILE__,__LINE__)
    this%elem_to_subset_id=0 
  end subroutine elem_to_subset_id_create

  !===============================================================================================
  subroutine elem_to_subset_id_free(this)
    implicit none
    class(elem_to_subset_id_t), intent(inout) :: this
    this%nelem= 0
    if (allocated(this%elem_to_subset_id)) call memfree (this%elem_to_subset_id,__FILE__,__LINE__)
  end subroutine elem_to_subset_id_free
  
  !===============================================================================================
  function elem_to_subset_id_get_subset_id ( this, ielem ) result ( elem_to_subset_id )
     class(elem_to_subset_id_t), intent(in) :: this
     integer(ip)               , intent(in) :: ielem  
     integer(ip)                            :: elem_to_subset_id
     elem_to_subset_id = this%elem_to_subset_id(ielem)
  end function elem_to_subset_id_get_subset_id

end module elem_to_subset_id_names
