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

!****************************************************************************************************

module abstract_reference_element_names
  use serial_names
# include "debug.i90"
  implicit none
  integer(ip), parameter :: DIM = number_space_dimensions
  private

  type object_tree_t
     private
     integer(ip) :: topology
     integer(ip), allocatable :: objects_array(:)     
     integer(ip), allocatable :: ijk_to_index(:)
     integer(ip) :: current = 1    
   contains
     procedure :: create_object_tree
     procedure, private :: fill_tree
     procedure :: create_children_iterator => object_tree_create_children_iterator
  end type object_tree_t

  type children_iterator_t
     private 
     type(object_tree_t), pointer :: object_tree
     integer(ip)                  :: parent
     integer(ip)                  :: component
     integer(ip)                  :: coordinate
   contains
     procedure :: create        => children_iterator_create     
     procedure :: current       => children_iterator_current
     procedure :: init          => children_iterator_init
     procedure :: next          => children_iterator_next
     procedure :: has_finished  => children_iterator_has_finished
!     procedure :: free          => children_iterator_free
     procedure, private :: is_admissible   
  end type children_iterator_t

  ! Types
  public :: object_tree_t, children_iterator_t

contains

  subroutine create_object_tree( this, topology )
    implicit none
    class(object_tree_t), intent(inout) :: this
    integer(ip)        , intent(in)    :: topology
    integer(ip) :: c, i

    this%topology = topology
    ! Compute number of object (all dimensions) for the hypercube
    ! c = (2**dim-1)*2**dim max id object
    c = 0
    do i = 0,DIM
       c = c + 2**(DIM-i)*(get_binomial_coefficient(DIM,i))
    end do
    ! Pre-allocate this%objects_array with c (exact for hypercube)
    call memalloc ( c, this%objects_array, __FILE__, __LINE__ )
    ! Pre-allocate the ijk_to_index. Maximum ijk_to_index is [111000] for dim = 3
    ! = ISHFT(ISHFT(1_ip,dim)-1,dim) = (2**dim-1)*2**dim
    call memalloc ( ISHFT(ISHFT(1_ip,DIM)-1,DIM), this%ijk_to_index, __FILE__, __LINE__  )
    this%ijk_to_index = 0
    ! Initialize current object 
    this%current = 0   
    ! Call the recursive fill_tree procedure starting with the volume object
    call this%fill_tree( ISHFT(ISHFT(1_ip,DIM)-1, DIM ) ) ! Root object (volume)
    ! Re-allocate the objects_array to the exact number of objects
    call memrealloc ( this%current, this%objects_array, __FILE__, __LINE__ )
    ! Sort objects based on its bit-based index
    !call sort( this%current, this%objects_array )
    ! Create tje bit-index to consecutive index array (0 for bit-indexes wo/ associated object)
    do i =1,this%current
       this%ijk_to_index( this%objects_array(i) ) = i
    end do
    ! Put the array pointing to volume object
    this%current = 1

  end subroutine create_object_tree

!=================================================================================================
! BNM(A,B)=A!/((A-B)!B!) computes the binomial coefficient of (A,B), A>B
integer (ip) function get_binomial_coefficient(a,b)
  implicit none
  integer(ip), intent(in)    :: a,b
  if (a >= b) then
     get_binomial_coefficient = int(get_factorial(a)/(get_factorial(b)*get_factorial(a-b)))
  else
     write(*,*) 'ERROR: no binomial coef for b > a'
     check(.false.)
  end if
end function get_binomial_coefficient

!==================================================================================================
! FC(i)=i! computes the factorial of i
integer(ip) function get_factorial(i)
  implicit none
  integer(ip), intent(in)    :: i
  integer(ip) :: k
  get_factorial = 1
  do k=2,i
     get_factorial = get_factorial*k
  end do
end function get_factorial

  recursive subroutine fill_tree( this, root )
    implicit none
    class(object_tree_t), intent(inout) :: this
    integer(ip)        , intent(in)    :: root
    type(children_iterator_t) :: children_iterator
    integer(ip)               :: children

    ! If the object not already inserted (use ijk_to_index as touch table)
    if ( this%ijk_to_index(root) /= -1 ) then 
       ! Increase one position current object in object_tree
       this%current = this%current + 1
       ! Put root as new object in object_tree
       this%objects_array(this%current) = root
       ! Mark root as touched (already inserted)
       this%ijk_to_index(root) = 1
       ! Create iterator over children of root object
       children_iterator = this%create_children_iterator(root)
       ! Loop over children
       do while (.not. children_iterator%has_finished() )
          children = children_iterator%current()
          ! Put childrens of child (recursive call)
          call this%fill_tree( children )
          ! Move to the next child
          call children_iterator%next()
       end do
    end if

  end subroutine fill_tree

  function object_tree_create_children_iterator(this, parent)
    implicit none
    class(object_tree_t), intent(in) :: this
    integer(ip)         , intent(in) :: parent
    type(children_iterator_t) :: object_tree_create_children_iterator
    call object_tree_create_children_iterator%create(this, parent)
  end function object_tree_create_children_iterator

  subroutine children_iterator_create ( this, object_tree, parent )
    implicit none
    class(children_iterator_t)        , intent(inout) :: this
    type(object_tree_t)       , target, intent(in)    :: object_tree
    integer(ip)                       , intent(in)    :: parent
    !call this%free()
    this%object_tree => object_tree
    this%parent = parent
    call this%init()
  end subroutine children_iterator_create

  subroutine children_iterator_init ( this )
    implicit none
    class(children_iterator_t), intent(inout) :: this
    this%component = 0
    this%coordinate = 0
    if ( .not. this%is_admissible() ) then
       call this%next()
    end if
  end subroutine children_iterator_init

  recursive subroutine children_iterator_next ( this )
    implicit none
    class(children_iterator_t), intent(inout) :: this
    if ( this%has_finished() ) return 
    if ( this%coordinate == 1 ) then
       this%component = this%component + 1
       this%coordinate = 0
    else
       this%coordinate = 1
    end if
    if ( .not. this%is_admissible() ) then
       call this%next()
    end if
  end subroutine children_iterator_next

  function children_iterator_has_finished ( this )
    implicit none
    class(children_iterator_t), intent(in) :: this
    logical :: children_iterator_has_finished
    children_iterator_has_finished = ( this%component >= DIM )
  end function children_iterator_has_finished

  function children_iterator_current ( this )
    implicit none
    class(children_iterator_t), intent(in) :: this
    integer(ip) :: children_iterator_current
    assert ( .not. this%has_finished() )
    children_iterator_current = get_child( this%parent, this%component, this%coordinate )
  end function children_iterator_current

  function is_admissible( this )
    implicit none
    class(children_iterator_t), intent(in) :: this
    logical                          :: is_admissible
    is_admissible = .false.
    if ( IBITS( this%parent, DIM + this%component, 1 ) == 1) then
       if ( IBITS( this%object_tree%topology, this%component, 1 ) == 1 ) then
          is_admissible = .true.
          ! Be careful w/ i = 0
       else if ( this%coordinate == 0 .or. this%component == 0 & 
            & .or. IBITS( this%parent, DIM, DIM+this%component-1 ) == 0 ) then
          is_admissible = .true.
       end if
    end if
  end function is_admissible

  function get_child( root, i, j )
    integer(ip), intent(in) :: root, i, j
    integer(ip) :: get_child
    get_child = root
    ! Fix component i
    get_child = IBCLR( get_child, DIM + i )
    ! Put value j in this component 
    if ( j == 1 ) then
       get_child = IBSET( get_child, i )
    end if
  end function get_child
end module abstract_reference_element_names

program abstract_reference_element
  use serial_names
  use abstract_reference_element_names
  implicit none
  type(object_tree_t) :: object_tree
  integer(ip), parameter :: DIM = number_space_dimensions

  ! 1D
  ! line 0 = 1 ( 0 = 1 )

  ! 2D
  ! square   10 = 11 ( 2 = 3 )
  ! triangle 00 = 01 ( 0 = 1 )

  ! 3D
  ! tetrahedron 000 = 001 ( 0 = 1 )
  ! pyramid     010 = 011 ( 2 = 3 )
  ! prysm       100 = 101 ( 4 = 5 )  
  ! cube        110 = 111 ( 6 = 7 )

  !call object_tree%create( 6 )

  write(*,*) 'DIMENSION',DIM
end program abstract_reference_element
