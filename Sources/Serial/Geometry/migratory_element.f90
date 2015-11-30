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
module migratory_element_names
  use types_names
  use hash_table_names
  use element_id_names
  implicit none
  private
# include "debug.i90"

  type, abstract :: migratory_element_t
     private
     class(element_id_t), allocatable :: id
   contains
     procedure (size_interface)  , deferred :: size
     procedure (pack_interface)  , deferred :: pack
     procedure (unpack_interface), deferred :: unpack
     procedure (free_interface)  , deferred :: free
     procedure (assignment_interface), deferred :: assign
     procedure :: get_id    => get_migratory_element_id
     procedure :: create_id => create_migratory_element_id
     generic   :: assignment(=) => assign
  end type migratory_element_t

  abstract interface
     subroutine size_interface(this,n)
       import :: migratory_element_t, ip
       implicit none
       class(migratory_element_t), intent(in)  :: this
       integer(ip)           , intent(out) :: n
     end subroutine size_interface

     subroutine pack_interface(this,n,buffer)
       import :: migratory_element_t, ip, ieep
       implicit none
       class(migratory_element_t), intent(in)  :: this
       integer(ip)             , intent(in)  :: n
       integer(ieep)           , intent(out) :: buffer(n)
     end subroutine pack_interface

     subroutine unpack_interface(this,n,buffer)
       import :: migratory_element_t, ip, ieep
       implicit none
       class(migratory_element_t), intent(inout) :: this
       integer(ip)             , intent(in)    :: n
       integer(ieep)           , intent(in)    :: buffer(n)
     end subroutine unpack_interface

     subroutine free_interface(this)
       import :: migratory_element_t
       implicit none
       class(migratory_element_t), intent(inout) :: this
     end subroutine free_interface

     subroutine assignment_interface(this,that)
       import :: migratory_element_t
       implicit none
       class(migratory_element_t), intent(inout) :: this
       class(migratory_element_t), intent(in)    :: that
     end subroutine assignment_interface

  end interface

  type, abstract :: migratory_element_pointer_t
     class(migratory_element_t), pointer :: p
  end type migratory_element_pointer_t

  public :: migratory_element_t, migratory_element_pointer_t

  !=============================================================================
  ! Abstract set and iterator
  !=============================================================================
  type, abstract :: migratory_element_iterator_t
     private
   contains
     procedure(next_interface)    , deferred :: next
     procedure(has_next_interface), deferred :: has_next
  end type migratory_element_iterator_t
  abstract interface
     function next_interface (this) result(p)
       import :: migratory_element_iterator_t, migratory_element_t
       implicit none
       class(migratory_element_iterator_t), intent(inout)  :: this
       class(migratory_element_t)         , pointer        :: p
     end function next_interface
     function has_next_interface(this) result(res)
       import :: migratory_element_iterator_t
       implicit none
       class(migratory_element_iterator_t), intent(inout) :: this
       logical :: res
     end function has_next_interface	
  end interface

  type, abstract :: migratory_element_set_t
   contains
     procedure(create_set_interface)     , deferred :: create
     procedure(free_set_interface)       , deferred :: free
     procedure(create_iterator_interface), deferred :: create_iterator
     procedure(free_iterator_interface)  , deferred :: free_iterator
  end type migratory_element_set_t
  ! Abstract interfaces
  abstract interface
     subroutine create_set_interface(this,size,mold)
       import :: migratory_element_set_t, ip, migratory_element_t
       implicit none
       class(migratory_element_set_t), intent(inout) :: this
       integer(ip)                  , intent(in)    :: size
       class(migratory_element_t)    , intent(in)    :: mold
     end subroutine create_set_interface
     subroutine free_set_interface(this)
       import :: migratory_element_set_t
       implicit none
       class(migratory_element_set_t), intent(inout) :: this
     end subroutine free_set_interface
     subroutine create_iterator_interface(this,iterator)
       import :: migratory_element_set_t, migratory_element_iterator_t
       implicit none
       class(migratory_element_set_t), target, intent(in)  :: this
       class(migratory_element_iterator_t)   , intent(out) :: iterator
     end subroutine create_iterator_interface
     subroutine free_iterator_interface(this,iterator)
       import :: migratory_element_set_t, migratory_element_iterator_t
       implicit none
       class(migratory_element_set_t)     , intent(in)    :: this
       class(migratory_element_iterator_t), intent(inout) :: iterator
     end subroutine free_iterator_interface
  end interface

! #define  template_element_t          migratory_element_t
! #define  template_element_set_t      migratory_element_set_t
! #define  template_element_iterator_t migratory_element_iterator_t
! #include "element_set.i90"
   public :: migratory_element_set_t, migratory_element_iterator_t

  !=============================================================================
  ! Plain set and iterator
  !=============================================================================

  type, extends(migratory_element_iterator_t) :: plain_migratory_element_iterator_t
     private
     integer(ip)                                 :: ielem=0
     type(plain_migratory_element_set_t), pointer :: plain_migratory_element_set
   contains
     procedure :: next     => plain_iterator_next
     procedure :: has_next => plain_iterator_has_next
  end type plain_migratory_element_iterator_t

  type, extends(migratory_element_set_t) :: plain_migratory_element_set_t
     private
     integer(ip)                            :: num_elements= -1  ! Number of elements in the set
     class(migratory_element_t), allocatable :: elements(:)       ! array of migratory_elements
   contains
     procedure :: create => create_plain_set
     procedure :: free   => free_plain_set
     procedure :: create_iterator => create_plain_iterator
     procedure :: free_iterator   => free_plain_iterator
  end type plain_migratory_element_set_t

! #define  plain_template_element_set_t      plain_migratory_element_set_t
! #define  plain_template_element_iterator_t plain_migratory_element_iterator_t
! #include "plain_element_set_header.i90"
   public :: plain_migratory_element_set_t, plain_migratory_element_iterator_t

  !=============================================================================
  ! Hash based set and iterator
  !=============================================================================

  type, extends(migratory_element_iterator_t) :: hash_migratory_element_iterator_t
     private
     class(element_id_t), allocatable  :: current_key
     class(element_id_t), allocatable  :: next_key
     type(hash_migratory_element_set_t), pointer :: hash_migratory_element_set
   contains
     procedure :: next     => hash_migratory_element_iterator_next
     procedure :: has_next => hash_migratory_element_iterator_has_next
  end type hash_migratory_element_iterator_t

  type hash_node
     private
     class(hash_node), pointer :: next_node => null()
     class(hash_node), pointer :: prev_node => null()
     !class(element_id_t)      , allocatable :: key
     class(migratory_element_t), allocatable :: element
   contains
     procedure :: put  => put_hash_node
     procedure :: get  => get_hash_node
     procedure :: get_key  => get_key_hash_node
     procedure :: get_next_key => get_next_key_hash_node
     procedure :: del  => del_hash_node
     procedure :: free => free_hash_node
  end type hash_node

  type, extends(migratory_element_set_t) :: hash_migratory_element_set_t
     private
     type(hash_node), dimension(:), allocatable :: vec
     integer(ip)                                :: vec_len = 0
     logical                                    :: is_init = .false.
   contains
     procedure :: create          => create_hash_migratory_element_set_t
     procedure :: free            => free_hash_migratory_element_set_t
     procedure :: create_iterator => create_hash_migratory_element_iterator
     procedure :: free_iterator   => free_hash_migratory_element_iterator
     procedure :: get             => get_hash_migratory_element_set_t
     procedure :: get_next_key    => get_next_key_hash_migratory_element_set_t
     procedure :: put             => put_hash_migratory_element_set_t
     procedure :: del             => del_hash_migratory_element_set_t
  end type hash_migratory_element_set_t

! #define  hash_template_element_set_t      hash_migratory_element_set_t
! #define  hash_template_element_iterator_t hash_migratory_element_iterator_t
! #include "hash_element_set_header.i90"
   public :: hash_migratory_element_set_t, hash_migratory_element_iterator_t

contains
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !=============================================================================
  function get_migratory_element_id(this) result(id)
    implicit none
    class(migratory_element_t), target, intent(in) :: this
    class(element_id_t)       , pointer    :: id
    id => this%id
  end function get_migratory_element_id
  subroutine create_migratory_element_id(this,id)
    implicit none
    class(migratory_element_t), intent(inout) :: this
    class(element_id_t)                       :: id
    if(allocated(this%id)) deallocate(this%id)
    allocate(this%id, mold=id)
  end subroutine create_migratory_element_id
  !=============================================================================
  !=============================================================================
  ! Plain set and iterator
  !=============================================================================
  !=============================================================================
  subroutine create_plain_set(this,size,mold)
    implicit none
    class(plain_migratory_element_set_t), intent(inout) :: this
    integer(ip)                         , intent(in)    :: size
    class(migratory_element_t)          , intent(in)    :: mold
    integer(ip) :: istat
    ! Create set
    this%num_elements = size
    allocate(this%elements(this%num_elements), mold=mold, stat=istat)
    check(istat==0)
  end subroutine create_plain_set

  subroutine free_plain_set(this)
    implicit none
    class(plain_migratory_element_set_t), intent(inout)  :: this
    integer(ip) :: istat
    this%num_elements = -1
    deallocate(this%elements, stat=istat)
    check(istat==0)
  end subroutine free_plain_set

  !=============================================================================
  subroutine create_plain_iterator(this,iterator)
    implicit none
    class(plain_migratory_element_set_t), target     , intent(in)  :: this
    class(migratory_element_iterator_t) , allocatable, intent(out) :: iterator
    integer(ip) :: istat
    allocate(plain_migratory_element_iterator_t :: iterator, stat=istat)
    check(istat==0)
    select type(iterator)
    class is(plain_migratory_element_iterator_t)
       iterator%plain_migratory_element_set => this
    end select
  end subroutine create_plain_iterator
  subroutine free_plain_iterator(this,iterator)
    implicit none
    class(plain_migratory_element_set_t)            , intent(in)    :: this
    class(migratory_element_iterator_t), allocatable, intent(inout) :: iterator
    ! No internal memory for this simple iterator.
    deallocate(iterator)
  end subroutine free_plain_iterator
  !=============================================================================
  function plain_iterator_next (this) result(p)
    implicit none
    class(plain_migratory_element_iterator_t), intent(inout) :: this
    class(migratory_element_t)               , pointer       :: p
    this%ielem = this%ielem + 1
    assert( this%ielem <= this%plain_migratory_element_set%num_elements )
    p => this%plain_migratory_element_set%elements(this%ielem)
  end function plain_iterator_next

  function plain_iterator_has_next(this) result(res)
    implicit none
    class(plain_migratory_element_iterator_t), intent(inout) :: this
    logical :: res
    res = ( this%ielem < this%plain_migratory_element_set%num_elements )
    if(.not.res) this%ielem=0 ! Reset iterator for the next loop
  end function plain_iterator_has_next

! #include "plain_element_set_body.i90"

  !=============================================================================
  !=============================================================================
  ! Hash based set and iterator
  !=============================================================================
  !=============================================================================
  ! node functions
  !
  recursive subroutine put_hash_node(list,element,stat)
    class(hash_node), target , intent(inout) :: list
    class(migratory_element_t), intent(in)   :: element
    integer(ip)              , intent(out)   :: stat
    class(hash_node)         , pointer       :: tmp => null()

    if( element%id%is_greater_than(list%element%id) ) then
       if(associated(list%next_node) ) then
          ! Keep going
          call put_hash_node(list%next_node,element,stat)
       else
          ! We are at the end of the list, allocate a new node
          allocate(list%next_node)
          list%next_node%prev_node => list
          allocate(list%element, mold=element)
          list%element = element
       end if
    else if( element%id%is_equal_to(list%element%id) ) then
       ! Already stored
       stat = was_stored
    else if( element%id%is_smaller_than(list%element%id) ) then
       ! Insert an element just after this one,
       tmp => list%next_node                                ! keep the reference to old next_node
       allocate(list%next_node)                             ! Allocate a new next_node
       allocate(list%next_node%element, mold=list%element)  ! WARNING: we need assignment working!
       list%next_node%element = list%element                ! New node fields (copy list keyval)
       list%next_node%prev_node => list
       list%next_node%next_node => tmp
       list%next_node%next_node%prev_node => list%next_node ! Regenerate link
       list%element = element                               ! Store keyval in list
       stat = now_stored
    end if
  end subroutine put_hash_node

  recursive subroutine del_hash_node(list,element,stat)
    class(hash_node)  , target, intent(inout) :: list
    class(migratory_element_t), intent(in)    :: element
    integer(ip)               , intent(out)   :: stat
    class(hash_node)          , pointer       :: tmp => null()

    if(element%id%is_equal_to(list%element%id)) then
       if(.not.associated(list%next_node)) then
          write(*,*) 'Error in migratory_element hash table'
          stat = error
          return
       end if
       ! WARNING: we need assignment working!
       list%element = list%next_node%element
       if(associated(list%next_node%next_node)) list%next_node%next_node%prev_node => list
       tmp => list%next_node%next_node
       call list%next_node%element%free()
       deallocate(list%next_node)
       list%next_node => tmp
       stat = deleted
    else if(associated(list%next_node)) then ! keep going
       call del_hash_node(list%next_node,element,stat)
    else
       stat = key_not_found
    end if

  end subroutine del_hash_node

  recursive subroutine get_hash_node(list,id,element,stat)
    class(hash_node)     , target     , intent(in)    :: list
    class(element_id_t)               , intent(in)    :: id
    class(migratory_element_t), pointer, intent(out)  :: element
    integer(ip)                       , intent(inout) :: stat
    if (id%is_equal_to(list%element%id)) then
       element => list%element
       stat = key_found
    else if(associated(list%next_node)) then ! keep going
       call get_hash_node(list%next_node,id,element,stat)
    else
       stat = key_not_found
    end if
  end subroutine get_hash_node

  recursive subroutine free_hash_node(list)
    implicit none
    class(hash_node), intent(inout) :: list
    if (associated(list%next_node)) then
       call free_hash_node(list%next_node)
       deallocate(list%next_node)
    end if
    call list%element%free()
    list%next_node => null()
  end subroutine free_hash_node

  ! Get the id of the next element in the list
  recursive subroutine get_next_key_hash_node(list,id,next_id,stat)
    class(hash_node)                  , intent(in)    :: list
    class(element_id_t)               , intent(in)    :: id
    class(element_id_t)               , intent(out)   :: next_id
    integer(ip)                       , intent(inout) :: stat
    if (id%is_equal_to(list%element%id)) then
       if(associated(list%next_node)) then
          next_id = list%next_node%element%id
          stat = key_found
       else
          stat = key_not_found
       end if
    else if(associated(list%next_node)) then
       call get_next_key_hash_node(list%next_node,id,next_id,stat)
    else
       stat = key_not_found
    end if
  end subroutine get_next_key_hash_node

  ! To get the key of the first element in the list
  subroutine  get_key_hash_node(list,id,stat)
    class(hash_node)                  , intent(in)    :: list
    class(element_id_t)               , intent(out)   :: id
    integer(ip)                       , intent(inout) :: stat
    if(associated(list%next_node)) then
       id = list%element%id
       stat = key_found
    else
       stat = key_not_found
    end if
  end subroutine get_key_hash_node

  !====================================================================
  ! tbl functions
  !
  subroutine create_hash_migratory_element_set_t(this,size,mold)
    class(hash_migratory_element_set_t), intent(inout) :: this
    integer(ip)                        , intent(in)    :: size
    class(migratory_element_t)         , intent(in)    :: mold
    this%vec_len = size
    if (allocated(this%vec)) deallocate(this%vec)
    allocate(this%vec(this%vec_len))
  end subroutine create_hash_migratory_element_set_t

  subroutine put_hash_migratory_element_set_t(tbl,element,stat)
    class(hash_migratory_element_set_t) , intent(inout) :: tbl
    class(migratory_element_t)          , intent(in)    :: element
    integer(ip)                         , intent(out)   :: stat
    integer(ip)                        :: hash
    hash = mod(element%id%to_int(),tbl%vec_len)
    assert ((hash>=1).and.(hash<=tbl%vec_len))
    call tbl%vec(hash)%put(element,stat)
  end subroutine put_hash_migratory_element_set_t

  subroutine get_hash_migratory_element_set_t(tbl,id,element,stat)
    class(hash_migratory_element_set_t), intent(in)    :: tbl
    class(element_id_t)                , intent(in)    :: id
    class(migratory_element_t), pointer, intent(out)   :: element
    integer(ip)                        , intent(out)   :: stat
    integer(ip)                                        :: hash
    hash = mod(id%to_int(),tbl%vec_len)
    assert ( (hash>=1).and.(hash<=tbl%vec_len))
    call tbl%vec(hash)%get(id,element,stat)
  end subroutine get_hash_migratory_element_set_t

  subroutine get_next_key_hash_migratory_element_set_t(tbl,key,next,stat)
    class(hash_migratory_element_set_t)  , intent(in) :: tbl
    class(element_id_t)                 , intent(in)  :: key
    class(element_id_t)                 , intent(out) :: next
    integer(ip)                         , intent(out) :: stat
    integer(ip)                                       :: hash
    hash = mod(key%to_int(),tbl%vec_len)
    assert ( (hash>=1).and.(hash<=tbl%vec_len))
    call tbl%vec(hash)%get_next_key(key,next,stat)
    do while(hash<=tbl%vec_len .and. stat == key_not_found )
       hash=hash+1
       call tbl%vec(hash)%get_key(next,stat)
    end do
  end subroutine get_next_key_hash_migratory_element_set_t

  subroutine del_hash_migratory_element_set_t(tbl,element,stat)
    class(hash_migratory_element_set_t), intent(inout) :: tbl
    class(migratory_element_t)         , intent(inout) :: element
    integer(ip)                        , intent(out)   :: stat
    integer(ip)                                        :: hash
    hash = mod(element%id%to_int(),tbl%vec_len)
    assert ( (hash>=1).and.(hash<=tbl%vec_len))
    call tbl%vec(hash)%del(element,stat)
  end subroutine del_hash_migratory_element_set_t

  subroutine free_hash_migratory_element_set_t(this)
    class(hash_migratory_element_set_t), intent(inout) :: this    
    integer(ip) :: i
    if (allocated(this%vec)) then
       do i=1,this%vec_len
          call this%vec(i)%free()
       end do
       deallocate(this%vec)
    end if
    this%is_init = .false.
  end subroutine free_hash_migratory_element_set_t

  !=============================================================================
  ! Simplest iterator
  subroutine create_hash_migratory_element_iterator(this,iterator)
    implicit none
    class(hash_migratory_element_set_t), target      , intent(in)  :: this
    class(migratory_element_iterator_t) , allocatable, intent(out) :: iterator
    integer(ip) :: istat
    allocate(hash_migratory_element_iterator_t :: iterator, stat=istat)
    check(istat==0)
    select type(iterator)
    class is(hash_migratory_element_iterator_t)
       iterator%hash_migratory_element_set => this
    end select
  end subroutine create_hash_migratory_element_iterator

  !=============================================================================
  subroutine free_hash_migratory_element_iterator(this,iterator)
    implicit none
    class(hash_migratory_element_set_t)            , intent(in)    :: this
    class(migratory_element_iterator_t), allocatable, intent(inout) :: iterator
    ! No internal memory for this simple iterator.
    deallocate(iterator)
  end subroutine free_hash_migratory_element_iterator

  !=============================================================================
  function hash_migratory_element_iterator_next (this) result(p)
    implicit none
    class(hash_migratory_element_iterator_t), intent(inout) :: this
    class(migratory_element_t)              , pointer       :: p
    integer(ip)  :: stat
    call this%hash_migratory_element_set%get( this%next_key, p, stat)
    assert(stat==key_found)
    this%current_key = this%next_key
  end function hash_migratory_element_iterator_next

  !=============================================================================
  function hash_migratory_element_iterator_has_next(this) result(res)
    implicit none
    class(hash_migratory_element_iterator_t), intent(inout) :: this
    logical :: res
    integer(ip)  :: stat
    call this%hash_migratory_element_set%get_next_key( this%current_key, this%next_key, stat)
    res = (stat==key_found)
    !if(.not.res) ! Reset iterator for the next loop to be done using default element_id
  end function hash_migratory_element_iterator_has_next
  
  ! #include "hash_element_set_body.i90"

end module migratory_element_names
