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
module JP_element_topology_names
  use types_names
  use memor_names
  use fe_space_types_names
  use element_id_names
  use migratory_element_names
  !use hash_table_names  
  implicit none
  private
# include "debug.i90"

  type, extends(migratory_element_t) :: JP_element_topology_t
     integer(ip)               :: num_vefs = -1    ! Number of vefs
     integer(ip), allocatable  :: vefs(:)          ! List of Local IDs of the vefs (vertices, edges, faces) that make up this element
     type(reference_element_t), pointer :: geo_reference_element => NULL() ! Topological info of the geometry (SBmod)
     real(rp), allocatable     :: coordinates(:,:)
     integer(ip)               :: order
   contains
     procedure :: size   => element_topology_size
     procedure :: pack   => element_topology_pack
     procedure :: unpack => element_topology_unpack
     procedure :: create => element_topology_create
     procedure :: free   => element_topology_free
     procedure :: assign => element_topology_assign
  end type JP_element_topology_t


  ! Types
  public :: JP_element_topology_t
  
  ! Functions
  public :: downcast_to_element_topology,local_id_from_vertices

contains

  subroutine element_topology_free(this)
    implicit none
    class(JP_element_topology_t), intent(inout) :: this
    if (allocated(this%vefs)) then
       call memfree(this%vefs, __FILE__, __LINE__)
    end if
    if (allocated(this%coordinates)) then
       call memfree(this%coordinates, __FILE__, __LINE__)
    end if
    this%num_vefs = -1
    nullify( this%geo_reference_element )
  end subroutine element_topology_free
 !=============================================================================
  subroutine element_topology_create(this)
    implicit none
    class(JP_element_topology_t), intent(inout) :: this
    !class(element_id_t)         , intent(in)    :: id_mold
    !call this%create_id(id_mold)
    assert(.not. allocated(this%vefs))
    this%num_vefs = -1
 end subroutine element_topology_create

 !=============================================================================
  function downcast_to_element_topology(parent) result(this)
    implicit none
    class(migratory_element_t), pointer, intent(in) :: parent
    class(JP_element_topology_t) , pointer             :: this
    select type(parent)
    class is(JP_element_topology_t)
       this => parent
    class default
       write(*,*) 'Cannot downcast to element_topology'
       check(.false.)
    end select
  end function downcast_to_element_topology

 !=============================================================================
  subroutine element_topology_assign(this, that)
    implicit none
    class(JP_element_topology_t), intent(inout) :: this
    class(migratory_element_t), intent(in)    :: that

    select type(that)
    class is(JP_element_topology_t)
       this = that
    class default
       write(*,*) 'Error calling JP_element_topology_t assignment'
       write(*,*) 'cannot assign object of another class'
       check(.false.)
    end select

  end subroutine element_topology_assign

 !=============================================================================
  subroutine local_id_from_vertices( e, nd, list, no, lid ) ! (SBmod)
    implicit none
    type(JP_element_topology_t), intent(in) :: e
    integer(ip), intent(in)  :: nd, list(:), no
    integer(ip), intent(out) :: lid
    ! Locals
    integer(ip)              :: first, last, io, iv, jv, ivl, c
    lid = -1

    do io = e%geo_reference_element%nvef_dim(nd), e%geo_reference_element%nvef_dim(nd+1)-1
       first =  e%geo_reference_element%crxob%p(io)
       last = e%geo_reference_element%crxob%p(io+1) -1
       if ( last - first + 1  == no ) then 
          do iv = first,last
             ivl = e%vefs(e%geo_reference_element%crxob%l(iv)) ! LID of vertices of the ef
             c = 0
             do jv = 1,no
                if ( ivl ==  list(jv) ) then
                   c  = 1 ! vertex in the external ef
                   exit
                end if
             end do
             if (c == 0) exit
          end do
          if (c == 1) then ! vef in the external element
             lid = e%vefs(io)
             exit
          end if
       end if
    end do
  end subroutine local_id_from_vertices

 !=============================================================================
 ! To be eliminated
 !=============================================================================
  subroutine element_topology_size (this, n)
    implicit none
    class(JP_element_topology_t), intent(in)  :: this
    integer(ip)            , intent(out) :: n
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip, size_of_rp

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_rp   = size(transfer(1.0_rp ,mold))
    n = 2*size_of_ip + size(this%coordinates)*size_of_rp

  end subroutine element_topology_size

  subroutine element_topology_pack (this, n, buffer)
    implicit none
    class(JP_element_topology_t), intent(in)  :: this
    integer(ip)              , intent(in)   :: n
    integer(ieep)            , intent(out)  :: buffer(n)
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip,size_of_rp
    integer(ip)   :: start, end

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_rp   = size(transfer(1.0_rp ,mold))

    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(size(this%coordinates,1),mold)

    start = end + 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(this%num_vefs,mold)

    !start = end + 1
    !end   = start + size_of_rp*size(this%coordinates) - 1
    !buffer(start:end) = transfer(this%coordinates,mold)

  end subroutine element_topology_pack

  subroutine element_topology_unpack(this, n, buffer)
    implicit none
    class(JP_element_topology_t), intent(inout)  :: this
    integer(ip)              , intent(in)     :: n
    integer(ieep)            , intent(in)     :: buffer(n)

    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip,size_of_rp
    integer(ip)   :: start, end
    integer(ip)   :: num_dims
    
    size_of_ip = size(transfer(1_ip ,mold))
    size_of_rp = size(transfer(1_rp,mold))

    start = 1
    end   = start + size_of_ip -1
    this%num_vefs  = transfer(buffer(start:end), this%num_vefs)

    start = end + 1
    end   = start + size_of_ip - 1
    num_dims  = transfer(buffer(start:end), num_dims)

   ! call memalloc( num_dims, this%num_vefs, this%coordinates, __FILE__, __LINE__ )
   ! start = end + 1
   ! end   = start + size_of_rp*size(this%coordinates) - 1
   ! this%coordinates  = transfer(buffer(start:end), this%coordinates)
    
  end subroutine element_topology_unpack

end module JP_element_topology_names
