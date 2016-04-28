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
module par_element_topology_names
  ! Serial modules
  use types_names
  use memor_names
  use hash_table_names
  use sort_names
  use migratory_element_names
  use cell_names
  !use hash_table_names

  ! Parallel modules
  !use par_context_names
  !use par_environment_names

  implicit none
# include "debug.i90"
  private

  type, extends(cell_t) :: par_element_topology_t
     !integer(ip)  :: interface  = -1           ! The boundary number ieboun (if this element is a interface element)
     integer(ip)  :: interface  = -1
     integer(ip)  :: mypart     = -1           ! To which part this element is mapped to ?
     integer(igp) :: globalID   = -1           ! Global ID of this element
                                               ! Local ID is given by the element_set
     integer(igp), allocatable :: vefs_GIDs(:) ! List of the GIDs of the vefs that make up this element
   contains
     procedure :: size     => par_element_topology_size
     procedure :: pack     => par_element_topology_pack
     procedure :: unpack   => par_element_topology_unpack
     procedure :: free     => par_element_topology_free
     procedure :: assign   => par_element_topology_assignment
     procedure :: add_data => par_element_topology_fill
     !procedure :: is_interface => par_element_topology_is_interface, par_element_topology_set_interface
     !procedure :: set_interface => par_element_topology_set_interface
  end type par_element_topology_t

  ! Types
  public :: par_element_topology_t

  public :: downcast_to_par_element_topology

  !=============================================================================
  ! Abstract set and iterator
  !=============================================================================
#define  abstract_element_t                        cell_t
#define  template_element_t                        par_element_topology_t
#define  abstract_template_element_set_t           par_element_topology_set_t
#define  abstract_template_element_iterator_t      par_element_topology_iterator_t
#include "abstract_element_set.i90"

  !=============================================================================
  ! Static set and iterator
  !=============================================================================
#define  static_template_element_set_t             static_par_element_topology_set_t
#define  static_template_element_iterator_t        static_par_element_topology_iterator_t
#include "static_element_set_header.i90"
   public :: static_cell_set_t, static_cell_iterator_t
 

contains

  !=============================================================================
  ! Static set and iterator
  !=============================================================================
#include "static_element_set_body.i90"

  !=============================================================================
  subroutine par_element_topology_create(this)
    implicit none
    class(par_element_topology_t), intent(inout) :: this

    call this%cell_t%create()
    assert(.not.allocated(this%vefs_GIDs))
    this%interface   = .false.
    this%mypart      = -1
    this%globalID    = -1
  end subroutine par_element_topology_create

  !=============================================================================
  subroutine par_element_topology_free(this)
    implicit none
    class(par_element_topology_t), intent(inout) :: this

    call this%cell_t%free()
    if (allocated(this%vefs_GIDs)) then
       call memfree(this%vefs_GIDs, __FILE__, __LINE__)
    end if
    this%interface = .false.
    this%mypart = -1 
  end subroutine par_element_topology_free
  
  !=============================================================================
  ! This function assumes that the data in element_topology has already been
  ! filled calling this%JP_element_topology_t%fill so it can be used separately.
  ! 
  subroutine par_element_topology_fill( this, globalID, mypart, vefs_GIDs)
    implicit none
    class(par_element_topology_t), intent(inout) :: this
    integer(igp), intent(in) :: globalID
    integer(ip) , intent(in) :: mypart
    integer(igp), intent(in) :: vefs_GIDs(:)

    call memalloc( this%num_vefs, this%vefs_GIDs, __FILE__, __LINE__)
    this%vefs_GIDs  = vefs_GIDs(1:this%num_vefs)
    this%mypart    = mypart
    this%globalID  = globalID

  end subroutine par_element_topology_fill

  !=============================================================================
  subroutine par_element_topology_assignment(this,that)
    implicit none
    class(par_element_topology_t)   , intent(inout) :: this
    class(migratory_element_t), intent(in)    :: that
    select type(that)
    class is(par_element_topology_t)
       this=that
    class default
       write(*,*) 'Error calling par_element_topology_t assignment'
       write(*,*) 'cannot assign object of another class'
       check(.false.)
    end select
  end subroutine par_element_topology_assignment

 !=============================================================================
  function downcast_to_par_element_topology(parent) result(this)
    implicit none
    class(migratory_element_t)    , pointer, intent(in) :: parent
    class(par_element_topology_t) , pointer             :: this
    select type(parent)
    class is(par_element_topology_t)
       this => parent
    class default
       write(*,*) 'Cannot downcast to element_topology'
       check(.false.)
    end select
  end function downcast_to_par_element_topology

  !=============================================================================
  subroutine par_element_topology_size (this, n)
    implicit none
    class(par_element_topology_t), intent(in)  :: this
    integer(ip)            , intent(out) :: n
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    n = size_of_ip*3 + size_of_igp*(this%num_vefs+1)

  end subroutine par_element_topology_size

  !=============================================================================
  subroutine par_element_topology_pack (this, n, buffer)
    implicit none
    class(par_element_topology_t), intent(in)  :: this
    integer(ip)            , intent(in)   :: n
    integer(ieep)            , intent(out)  :: buffer(n)
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    integer(ip) :: start, end

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(this%mypart,mold)

    start = end + 1
    end   = start + size_of_ip - 1
    buffer(start:end) = transfer(this%interface,mold)

    start = end + 1
    end   = start + size_of_igp - 1 
    buffer(start:end) = transfer(this%globalID,mold)

    start = end + 1
    end   = start + size_of_ip - 1
    buffer(start:end) = transfer(this%num_vefs,mold)

    start = end + 1
    end   = start + this%num_vefs*size_of_igp - 1
    buffer(start:end) = transfer(this%vefs_GIDs,mold)

  end subroutine par_element_topology_pack

  !=============================================================================
  subroutine par_element_topology_unpack(this, n, buffer)
    implicit none
    class(par_element_topology_t), intent(inout) :: this
    integer(ip)            , intent(in)     :: n
    integer(ieep)            , intent(in)     :: buffer(n)

    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    integer(ip) :: start, end
    
    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))
    
    start = 1
    end   = start + size_of_ip -1
    this%mypart  = transfer(buffer(start:end), this%mypart)

    start = end + 1
    end   = start + size_of_ip - 1
    this%interface  = transfer(buffer(start:end), this%interface)

    start = end + 1
    end   = start + size_of_igp - 1 
    this%globalID = transfer(buffer(start:end), this%globalID)

    start = end + 1
    end   = start + size_of_ip - 1
    this%num_vefs = transfer(buffer(start:end), this%num_vefs)

    call memalloc( this%num_vefs, this%vefs_GIDs, __FILE__, __LINE__ )

    start = end + 1
    end   = start + this%num_vefs*size_of_igp - 1
    this%vefs_GIDs = transfer(buffer(start:end), this%vefs_GIDs)

  end subroutine par_element_topology_unpack

end module par_element_topology_names
