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
module cell_names
  use types_names
  use memor_names
  use list_types_names
  use hash_table_names
  !use element_id_names
  use migratory_element_names
  use reference_fe_names
  implicit none
  private
# include "debug.i90"

  ! Explicit control of element states is NOT implemented
  ! for this class (because it requires an integer).
  ! Conceptually we consider the following diagram:
  !
  ! --------------------------------------------------------
  ! Input State      | Action               | Output State 
  ! --------------------------------------------------------
  ! not_created      | create               | created 
  ! not_created      | free                 | not_created
  ! created          | free                 | not_created
  ! created          | fill                 | filled
  ! filled           | fill                 | filled
  ! filled           | free                 | not_created

  type, extends(migratory_element_t) :: cell_t
     integer(ip)               :: num_vefs = -1    ! Number of vefs
     integer(ip), allocatable  :: vefs(:)          ! List of Local IDs of the vefs (vertices, edges, faces) that make up this element
     real(rp), allocatable     :: coordinates(:,:)
     integer(ip)               :: order
     class(reference_fe_t), pointer :: geo_reference_element => NULL() ! Topological info of the geometry (SBmod)
   contains
     procedure :: create => cell_create
     procedure :: free   => cell_free
     procedure :: assign => cell_assign
     procedure :: fill   => cell_fill
     procedure :: print  => cell_print
  end type cell_t

  ! Types
  public :: cell_t
  
  ! Functions
  public :: downcast_to_cell, local_id_from_vertices

  !=============================================================================
  ! Abstract set and iterator, just needed to define the static, not used as such.
  !=============================================================================
#define  template_element_t                        cell_t
#define  abstract_template_element_set_t           abstract_cell_set_t
#define  abstract_template_element_iterator_t      abstract_cell_iterator_t
#include "abstract_element_set.i90"

  !=============================================================================
  ! Static set and iterator
  !=============================================================================
#define  static_template_element_set_t             static_cell_set_t
#define  static_template_element_iterator_t        static_cell_iterator_t
#include "static_element_set_header.i90"
   public :: static_cell_set_t, static_cell_iterator_t

  !=============================================================================
  ! Hash based set and iterator
  !=============================================================================
! #define  key_type                         integer(ip)
! #define  hash_template_element_set_t      hash_cell_set_t
! #define  hash_template_element_iterator_t hash_cell_iterator_t
! #include "hash_element_set_header.i90"
!    public :: hash_cell_set_t, hash_cell_iterator_t
  
contains

 !=============================================================================
  subroutine cell_create(this)
    implicit none
    class(cell_t), intent(inout) :: this
    assert(.not. allocated(this%vefs))
    assert(.not. allocated(this%coordinates))
    this%num_vefs = -1
 end subroutine cell_create

 !=============================================================================
  subroutine cell_free(this)
    implicit none
    class(cell_t), intent(inout) :: this
    if (allocated(this%vefs)) then
       call memfree(this%vefs, __FILE__, __LINE__)
    end if
    if (allocated(this%coordinates)) then
       call memfree(this%coordinates, __FILE__, __LINE__)
    end if
    this%num_vefs = -1
    nullify( this%geo_reference_element )
  end subroutine cell_free

 !=============================================================================
  subroutine cell_fill(this, num_dims, num_vefs, num_nodes, vefs, coordinates)
    implicit none
    class(cell_t), intent(inout) :: this
    integer(ip), intent(in) :: num_vefs, num_dims, num_nodes
    integer(ip), intent(in) :: vefs(num_vefs)
    real(rp)   , intent(in) :: coordinates(:,:)

    this%num_vefs = num_vefs
    call memalloc( this%num_vefs, this%vefs, __FILE__, __LINE__)
    this%vefs = vefs
    call memalloc( num_dims, num_nodes, this%coordinates, __FILE__, __LINE__ )
    this%coordinates = coordinates

  end subroutine cell_fill

 !=============================================================================
  subroutine cell_assign(this, that)
    implicit none
    class(cell_t), intent(inout) :: this
    class(migratory_element_t), intent(in)    :: that

    select type(that)
    class is(cell_t)
       this = that
    class default
       write(*,*) 'Error calling JP_cell_t assignment'
       write(*,*) 'cannot assign object of another class'
       check(.false.)
    end select

  end subroutine cell_assign

 !=============================================================================
  subroutine cell_print(elem, lunou)
    implicit none
    class(cell_t), intent(inout) :: elem
    integer(ip)                 , intent(in)    :: lunou
    write (lunou,*) 'num_vefs:'   , elem%num_vefs
    write (lunou,*) 'vefs:'       , elem%vefs
    write (lunou,*) 'coordinates:', elem%coordinates
    write (lunou,*) 'order:'      , elem%order
    !call reference_element_write ( elem%geo_reference_element )
  end subroutine cell_print

 !=============================================================================
  function downcast_to_cell(parent) result(this)
    implicit none
    class(migratory_element_t), pointer, intent(in) :: parent
    class(cell_t) , pointer             :: this
    select type(parent)
    class is(cell_t)
       this => parent
    class default
       write(*,*) 'Cannot downcast to cell'
       check(.false.)
    end select
  end function downcast_to_cell

 !=============================================================================
  !subroutine local_id_from_vertices( e, nd, list, no, lid ) ! (SBmod)
  !  implicit none
  !  class(JP_cell_t), intent(in) :: e
  !  integer(ip), intent(in)  :: nd, list(:), no
  !  integer(ip), intent(out) :: lid
  !  ! Locals
  !  integer(ip)              :: first, last, io, iv, jv, ivl, c
  !  lid = -1

  !  do io = e%geo_reference_element%nvef_dim(nd), e%geo_reference_element%nvef_dim(nd+1)-1
  !     first =  e%geo_reference_element%crxob%p(io)
  !     last = e%geo_reference_element%crxob%p(io+1) -1
  !     if ( last - first + 1  == no ) then 
  !        do iv = first,last
  !           ivl = e%vefs(e%geo_reference_element%crxob%l(iv)) ! LID of vertices of the ef
  !           c = 0
  !           do jv = 1,no
  !              if ( ivl ==  list(jv) ) then
  !                 c  = 1 ! vertex in the external ef
  !                 exit
  !              end if
  !           end do
  !           if (c == 0) exit
  !        end do
  !        if (c == 1) then ! vef in the external element
  !           lid = e%vefs(io)
  !           exit
  !        end if
  !     end if
  !  end do
  !end subroutine local_id_from_vertices

   subroutine local_id_from_vertices( e, nd, list, no, lid ) ! (SBmod)
    implicit none
    class(cell_t), intent(in) :: e
    integer(ip), intent(in)  :: nd, list(:), no
    integer(ip), intent(out) :: lid
    ! Locals
    integer(ip)              :: first, last, io, iv, jv, ivl, c
    type(list_t), pointer :: vertices_vef
    
    vertices_vef => e%geo_reference_element%get_vertices_vef()
    lid = -1
    do io = e%geo_reference_element%get_first_vef_id_of_dimension(nd-1), e%geo_reference_element%get_first_vef_id_of_dimension(nd)-1
       first =  vertices_vef%p(io)
       last = vertices_vef%p(io+1) -1
       if ( last - first + 1  == no ) then 
          do iv = first,last
             ivl = e%vefs(vertices_vef%l(iv)) ! LID of vertices of the ef
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
  ! subroutine element_topology_size (this, n)
  !   implicit none
  !   class(JP_element_topology_t), intent(in)  :: this
  !   integer(ip)            , intent(out) :: n
    
  !   ! Locals
  !   integer(ieep) :: mold(1)
  !   integer(ip)   :: size_of_ip, size_of_rp

  !   size_of_ip   = size(transfer(1_ip ,mold))
  !   size_of_rp   = size(transfer(1.0_rp ,mold))
  !   n = 2*size_of_ip + size(this%coordinates)*size_of_rp

  ! end subroutine element_topology_size

  ! subroutine element_topology_pack (this, n, buffer)
  !   implicit none
  !   class(JP_element_topology_t), intent(in)  :: this
  !   integer(ip)              , intent(in)   :: n
  !   integer(ieep)            , intent(out)  :: buffer(n)
    
  !   ! Locals
  !   integer(ieep) :: mold(1)
  !   integer(ip)   :: size_of_ip,size_of_rp
  !   integer(ip)   :: start, end

  !   size_of_ip   = size(transfer(1_ip ,mold))
  !   size_of_rp   = size(transfer(1.0_rp ,mold))

  !   start = 1
  !   end   = start + size_of_ip -1
  !   buffer(start:end) = transfer(size(this%coordinates,1),mold)

  !   start = end + 1
  !   end   = start + size_of_ip -1
  !   buffer(start:end) = transfer(this%num_vefs,mold)

  !   !start = end + 1
  !   !end   = start + size_of_rp*size(this%coordinates) - 1
  !   !buffer(start:end) = transfer(this%coordinates,mold)

  ! end subroutine element_topology_pack

  ! subroutine element_topology_unpack(this, n, buffer)
  !   implicit none
  !   class(JP_element_topology_t), intent(inout)  :: this
  !   integer(ip)              , intent(in)     :: n
  !   integer(ieep)            , intent(in)     :: buffer(n)

  !   ! Locals
  !   integer(ieep) :: mold(1)
  !   integer(ip)   :: size_of_ip,size_of_rp
  !   integer(ip)   :: start, end
  !   integer(ip)   :: num_dims
    
  !   size_of_ip = size(transfer(1_ip ,mold))
  !   size_of_rp = size(transfer(1_rp,mold))

  !   start = 1
  !   end   = start + size_of_ip -1
  !   this%num_vefs  = transfer(buffer(start:end), this%num_vefs)

  !   start = end + 1
  !   end   = start + size_of_ip - 1
  !   num_dims  = transfer(buffer(start:end), num_dims)

  !  ! call memalloc( num_dims, this%num_vefs, this%coordinates, __FILE__, __LINE__ )
  !  ! start = end + 1
  !  ! end   = start + size_of_rp*size(this%coordinates) - 1
  !  ! this%coordinates  = transfer(buffer(start:end), this%coordinates)
    
  ! end subroutine element_topology_unpack

  !=============================================================================
  !=============================================================================
  ! Static set and iterator
  !=============================================================================
  !=============================================================================
#include "static_element_set_body.i90"

  !=============================================================================
  !=============================================================================
  ! Hash based set and iterator
  !=============================================================================
  !=============================================================================
!#define convert_to_int(key) key
!#include "hash_element_set_body.i90"
  
  
end module cell_names
