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

!module abstract_reference_element_names
!  use serial_names
!  use sort_names
!# include "debug.i90"
!  implicit none
!  integer(ip), parameter :: DIM = 2
!  private



!  type geometry_tree_t
!     !private
!     integer(ip)              :: topology
!     integer(ip)              :: number_objects 
!     integer(ip), allocatable :: object_array(:)     
!     integer(ip), allocatable :: ijk_to_index(:)
!   contains
!     procedure          :: create                   => geometry_tree_create
!     procedure          :: create_children_iterator => geometry_tree_create_children_iterator
!     procedure, private :: fill_tree 
!  end type geometry_tree_t

!  type node_array_t
!     !private
!     type(geometry_tree_t), pointer :: object_tree
!     integer(ip)                  :: order
!     integer(ip)                  :: number_nodes
!     integer(ip), allocatable     :: node_array(:)
!     integer(ip), allocatable     :: ijk_to_index(:)
!   contains
!     procedure :: create               => node_array_create
!     procedure :: fill                 => node_array_fill
!     procedure :: print                => node_array_print
!     !procedure :: free          => node_array_free
!     procedure :: create_node_iterator => node_array_create_node_iterator
!  end type node_array_t

!  type children_iterator_t
!     private 
!     type(geometry_tree_t), pointer :: object_tree
!     integer(ip)                  :: parent
!     integer(ip)                  :: component
!     integer(ip)                  :: coordinate
!   contains
!     procedure :: create        => children_iterator_create     
!     procedure :: current       => children_iterator_current
!     procedure :: init          => children_iterator_init
!     procedure :: next          => children_iterator_next
!     procedure :: has_finished  => children_iterator_has_finished
!     !     procedure :: free          => children_iterator_free
!     procedure :: print         => children_iterator_print
!     procedure, private :: is_admissible => children_iterator_is_admissible   
!  end type children_iterator_t

!  type node_iterator_t
!     !private 
!     type(node_array_t), pointer :: node_array
!     logical                     :: own_boundary
!     integer(ip)                 :: object
!     integer(ip)                 :: topology
!     integer(ip)                 :: displacement(0:DIM-1)
!     integer(ip)                 :: coordinate(0:DIM-1)
!     logical                     :: overflow
!     integer(ip)                  :: max_value ! 0 or 1
!     integer(ip)                  :: min_value ! order or order-1
!   contains
!     procedure :: create        => node_iterator_create     
!     procedure :: current       => node_iterator_current
!     procedure :: init          => node_iterator_init
!     procedure :: next          => node_iterator_next
!     procedure :: has_finished  => node_iterator_has_finished
!     !procedure :: free          => node_iterator_free
!     procedure :: print         => node_iterator_print
!     procedure, private :: in_bound => node_iterator_in_bound 
!  end type node_iterator_t

!  ! Types
!  public :: geometry_tree_t, node_array_t, children_iterator_t, node_iterator_t

!contains

!  !=================================================================================================
!  ! This subroutine creates the object tree, where the geometrical objects are defined by 2*DIM bits.
!  ! The "first" (from the right) DIM bits provide the "anchor node", whereas the "second" set of DIM
!  ! bits provides the type of object. For every component, if we have 1 it means that the y-component
!  ! can freely change. When = 0 it is constrained. However, it does not mean that it is constant (it
!  ! is only for hyper-cubes). Think about the face of the tetrahedron that is not aligned with any 
!  ! axes.
!  subroutine geometry_tree_create( this, topology )
!    implicit none
!    class(geometry_tree_t), intent(inout) :: this
!    integer(ip)        , intent(in)    :: topology
!    integer(ip) :: c, i, istat
!    this%topology = topology
!    ! Compute number of object (all dimensions) for the hypercube
!    ! c = (2**dim-1)*2**dim max id object
!    c = 0
!    do i = 0,DIM
!       c = c + 2**(DIM-i)*(get_binomial_coefficient(DIM,i))
!    end do
!    ! Pre-allocate this%objects_array with c (exact for hypercube)
!    call memalloc ( c, this%object_array, __FILE__, __LINE__ )
!    ! Pre-allocate the ijk_to_index. Maximum ijk_to_index is [111000]+1 for dim = 3
!    ! = ISHFT(ISHFT(1_ip,dim)-1,dim) = (2**dim-1)*2**dim
!    call memalloc ( ISHFT(ISHFT(1_ip,DIM)-1,DIM)+1, this%ijk_to_index, __FILE__, __LINE__, lb1 = 0 )
!    this%ijk_to_index = 0
!    this%number_objects = 0
!    ! Call the recursive fill_tree procedure starting with the volume object
!    call this%fill_tree( ISHFT(ISHFT(1_ip,DIM)-1, DIM ) ) ! Root object (volume)
!    ! Re-allocate the objects_array to the exact number of objects
!    call memrealloc ( this%number_objects, this%object_array, __FILE__, __LINE__ )
!    ! Sort objects based on its bit-based index
!    call sort( this%number_objects, this%object_array )
!    ! Create the bit-index to consecutive index array (0 for bit-indexes wo/ associated object)
!    do i =1,this%number_objects
!       this%ijk_to_index( this%object_array(i) ) = i
!    end do
!  end subroutine geometry_tree_create

!  !=================================================================================================
!  ! get_binomial_coefficient(A,B)=A!/((A-B)!B!) computes the binomial coefficient of (A,B), A>B
!  integer (ip) function get_binomial_coefficient(a,b)
!    implicit none
!    integer(ip), intent(in)    :: a,b
!    if (a >= b) then
!       get_binomial_coefficient = int(get_factorial(a)/(get_factorial(b)*get_factorial(a-b)))
!    else
!       write(*,*) 'ERROR: no binomial coef for b > a'
!       check(.false.)
!    end if
!  end function get_binomial_coefficient

!  !==================================================================================================
!  ! get_factorial(i)=i! computes the factorial of i
!  integer(ip) function get_factorial(i)
!    implicit none
!    integer(ip), intent(in)    :: i
!    integer(ip) :: k
!    get_factorial = 1
!    do k=2,i
!       get_factorial = get_factorial*k
!    end do
!  end function get_factorial

!  !==================================================================================================
!  ! This recursive subroutine fills ALL the sub-objects (including the object itself). From dimension
!  ! DIM to 0.
!  recursive subroutine fill_tree( this, root )
!    implicit none
!    class(geometry_tree_t), intent(inout) :: this
!    integer(ip)        , intent(in)    :: root
!    type(children_iterator_t) :: children_iterator
!    integer(ip)               :: children  
!    ! If the object not already inserted (use ijk_to_index as touch table)
!    
!    if ( this%ijk_to_index(root) == 0 ) then 
!       ! Increase one position current object in object_tree
!       this%number_objects = this%number_objects + 1
!       write(*,*) 'NEW OBJECT IN POS', this%number_objects
!       ! Put root as new object in object_tree
!       this%object_array(this%number_objects) = root
!       write(*,*) 'OBJECT INDEX:'
!       write(*,'(B32)') root
!       ! Mark root as touched (already inserted)
!       this%ijk_to_index(root) = 1
!       ! Create iterator over children of root object
!       children_iterator = this%create_children_iterator(root)
!       ! Loop over children
!       do while (.not. children_iterator%has_finished() )
!          children = children_iterator%current()
!          !write(*,*) 'current children: '
!          !write(*,'(B32)') children
!          ! Put children of child (recursive call)
!          call this%fill_tree( children )
!          ! Move to the next child
!          call children_iterator%next()
!       end do
!    end if

!  end subroutine fill_tree

!  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  
!  !==================================================================================================
!  function geometry_tree_create_children_iterator(this, parent)
!    implicit none
!    class(geometry_tree_t), intent(in) :: this
!    integer(ip)         , intent(in) :: parent
!    type(children_iterator_t) :: geometry_tree_create_children_iterator
!    call geometry_tree_create_children_iterator%create(this, parent)
!  end function geometry_tree_create_children_iterator

!  subroutine children_iterator_create ( this, object_tree, parent )
!    implicit none
!    class(children_iterator_t)        , intent(inout) :: this
!    type(geometry_tree_t)       , target, intent(in)    :: object_tree
!    integer(ip)                       , intent(in)    :: parent
!    !call this%free()
!    this%object_tree => object_tree
!    this%parent = parent
!    call this%init()
!  end subroutine children_iterator_create

!  subroutine children_iterator_init ( this )
!    implicit none
!    class(children_iterator_t), intent(inout) :: this
!    this%component = 0
!    this%coordinate = 0
!    if ( .not. this%is_admissible() ) then
!       call this%next()
!    end if
!  end subroutine children_iterator_init

!  subroutine children_iterator_print( this )
!    implicit none
!    class(children_iterator_t), intent(inout) :: this
!    write(*,*) '***CHILDREN_ITERATOR ***'
!    write(*,*) 'parent: '
!    write(*,'(B32)') this%parent
!    write(*,*) 'component: ',this%component
!    write(*,*) 'coordinate: ',this%coordinate

!  end subroutine children_iterator_print

!  recursive subroutine children_iterator_next ( this )
!    implicit none
!    class(children_iterator_t), intent(inout) :: this
!    if ( this%has_finished() ) return 
!    if ( this%coordinate == 1 ) then
!       this%component = this%component + 1
!       this%coordinate = 0
!    else
!       this%coordinate = 1
!    end if
!    if ( .not. this%is_admissible() ) then
!       call this%next()
!    else 
!    end if
!  end subroutine children_iterator_next

!  function children_iterator_has_finished ( this )
!    implicit none
!    class(children_iterator_t), intent(in) :: this
!    logical :: children_iterator_has_finished
!    children_iterator_has_finished = ( this%component >= DIM )
!  end function children_iterator_has_finished

!  function children_iterator_current ( this )
!    implicit none
!    class(children_iterator_t), intent(in) :: this
!    integer(ip) :: children_iterator_current, j
!    assert ( .not. this%has_finished() )
!    children_iterator_current = this%parent
!    children_iterator_current = IBCLR( children_iterator_current, DIM + this%component )
!    if ( this%coordinate == 1 ) then
!       if ( IBITS( this%object_tree%topology, this%component, 1 )  == 0 ) then
!          do j = 0,this%component-1
!             children_iterator_current = IBCLR( children_iterator_current, j )
!          end do
!       end if
!       children_iterator_current = IBSET( children_iterator_current, this%component )
!    end if
!  end function children_iterator_current

!  function children_iterator_is_admissible( this )
!    implicit none
!    class(children_iterator_t), intent(inout) :: this
!    logical                          :: children_iterator_is_admissible
!    children_iterator_is_admissible = .false. 
!    if ( IBITS( this%parent, DIM + this%component, 1 ) == 1) then
!       if ( IBITS( this%object_tree%topology, this%component, 1 ) == 1 ) then
!          children_iterator_is_admissible = .true.
!       else if ( this%coordinate == 0 .or. this%component == 0 & 
!            & .or. IBITS( this%parent, DIM, this%component ) == 0 ) then  
!          children_iterator_is_admissible = .true.
!       end if
!    end if
!  end function children_iterator_is_admissible

!  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  function node_array_create_node_iterator( node_array, parent, own_boundary)
!    implicit none
!    class(node_array_t) , intent(in) :: node_array
!    integer(ip)         , intent(in) :: parent
!    logical             , intent(in) :: own_boundary 
!    type(node_iterator_t) :: node_array_create_node_iterator
!    call node_array_create_node_iterator%create( node_array, parent, own_boundary )
!  end function node_array_create_node_iterator

!  subroutine node_iterator_create ( this, node_array, object, own_boundary )
!    implicit none
!    class(node_iterator_t)            , intent(inout) :: this
!    type(node_array_t), target               , intent(in)    :: node_array
!    integer(ip)                       , intent(in)    :: object
!    logical                           , intent(in)    :: own_boundary
!    integer(ip) :: c, i
!    this%node_array => node_array
!    this%topology = 0
!    this%object = object
!    this%own_boundary = own_boundary
!    if ( own_boundary) then
!       this%max_value = node_array%order
!       this%min_value = 0
!    else
!       this%max_value = node_array%order-1
!       this%min_value = 1
!    end if
!    ! Restrict the topology of the object_tree to the current object
!    do i = 0, DIM-1
!       if ( IBITS( this%object, DIM+i, 1 ) == 1 ) then 
!          if ( IBITS( this%node_array%object_tree%topology, i, 1 ) == 1 ) this%topology = IBSET( this%topology, c )
!          c = c+1
!       end if
!    end do
!  end subroutine node_iterator_create

!  subroutine node_iterator_init ( this )
!    implicit none
!    class(node_iterator_t), intent(inout) :: this
!    integer(ip) :: c, i
!    this%displacement = 0
!    this%coordinate = 0
!    this%overflow = .true.
!    if ( object_dimension(this%object) > 0 .or. this%own_boundary ) this%overflow = .false.
!    c = 0
!    do i = 0, DIM-1
!       if ( IBITS( this%object, DIM+i, 1 ) == 0 ) then  ! if is fixed component in object
!          if ( IBITS( this%object, i, 1 ) == 0 ) then   ! if coordinate of anchor node == 0
!             this%coordinate(i) = 0                     ! fixed coordinate = min value
!          else
!             this%coordinate(i) = this%node_array%order ! fixed coordinate = max value   
!          end if
!       else 
!          this%displacement(c) = this%min_value
!          c = c+1
!       end if
!    end do
!    if ( object_dimension(this%object) > 0 ) then
!       this%displacement(0) = this%displacement(0)-1
!       call this%next()
!    end if
!  end subroutine node_iterator_init

!  function object_dimension( object )
!    integer(ip), intent(in) :: object
!    integer(ip) :: object_dimension, i
!    object_dimension = 0
!    do i = 0,DIM-1
!       object_dimension = object_dimension + IBITS( object, DIM+i, 1 )
!    end do
!  end function object_dimension

!  subroutine node_iterator_print( this )
!    implicit none
!    class(node_iterator_t), intent(inout) :: this
!    write(*,*) '***NODE_ITERATOR ***'
!    write(*,*) 'object: '
!    write(*,'(B32)') this%object
!    write(*,*) 'displacement: ',this%displacement
!  end subroutine node_iterator_print

!  subroutine node_iterator_next ( this )
!    implicit none
!    class(node_iterator_t), intent(inout) :: this
!    integer(ip) :: comp, end_comp, i
!    if ( this%has_finished() ) return 
!    comp = 0 
!    end_comp = object_dimension(this%object)-1
!    this%overflow = .true.
!    do while( comp <= end_comp)
!       if ( this%overflow ) then 
!          this%displacement(comp) = this%displacement(comp)+1
!          if ( this%in_bound(comp, end_comp) ) then
!             this%overflow = .false.
!             do i = comp-1,0,-1
!                if ( .not. this%in_bound( i, end_comp ) ) then
!                   this%overflow = .true.
!                end if
!             end do
!             exit
!          else 
!             this%displacement(comp) = this%min_value
!          end if
!       end if
!       comp = comp+1
!    end do
!  end subroutine node_iterator_next

!  function node_iterator_in_bound( this, comp, end_comp )
!    implicit none
!    class(node_iterator_t), intent(in) :: this
!    integer(ip)           , intent(in) :: comp, end_comp
!    logical :: node_iterator_in_bound
!    integer(ip) :: bound, i
!    node_iterator_in_bound = .true.
!    bound = this%max_value
!    do i = comp+1,end_comp
!       if ( IBITS( this%topology, i, 1 ) == 0 ) then
!          bound = bound - this%displacement(i)
!       end if
!    end do
!    if ( this%displacement(comp) > bound ) node_iterator_in_bound = .false.
!  end function node_iterator_in_bound

!  function node_iterator_has_finished ( this )
!    implicit none
!    class(node_iterator_t), intent(in) :: this
!    logical :: node_iterator_has_finished
!    node_iterator_has_finished = this%overflow
!  end function node_iterator_has_finished

!  function node_iterator_current ( this )
!    implicit none
!    class(node_iterator_t), intent(inout) :: this
!    integer(ip) :: node_iterator_current
!    integer(ip) :: coord(0:DIM-1), c, comp, j
!    assert ( .not. this%has_finished() )
!    ! First put free coordinates as such (fixed ones introduced when initializing)
!    c = 0
!    do comp =0,DIM-1
!       if ( IBITS ( this%object, comp+DIM, 1 ) == 1 ) then
!          this%coordinate(comp) = this%displacement(c)
!          c = c+1
!       end if
!    end do
!    ! Next, translate fixed coordinates if needed
!    do comp =0,DIM-1
!       if ( IBITS ( this%object, comp+DIM, 1 ) == 0 .and. IBITS( this%object, comp, 1 ) == 1 ) then
!          this%coordinate(comp) = this%node_array%order
!          do j = comp+1,DIM-1
!             if ( IBITS ( this%node_array%object_tree%topology, j, 1 ) == 0 .and. IBITS ( this%object, j+DIM, 1 ) == 1 ) then
!                this%coordinate(comp) = this%coordinate(comp) - this%coordinate(j)
!             end if
!          end do
!       end if
!    end do
!    node_iterator_current = ijk_to_index( this%coordinate, this%node_array%order )
!  end function node_iterator_current

!  function ijk_to_index( n, order )
!    integer(ip), intent(in) :: n(0:DIM-1), order
!    integer(ip) :: ijk_to_index, i
!    ijk_to_index = 0
!    do i = 0,DIM-1
!       ijk_to_index = ijk_to_index + n(i)*((order+1)**i)
!    end do
!  end function ijk_to_index

!  subroutine node_array_create ( this, object_tree, order )
!    implicit none
!    class(node_array_t)        , intent(inout) :: this
!    type(geometry_tree_t), target, intent(in)    :: object_tree
!    integer(ip)                , intent(in)    :: order
!    this%object_tree => object_tree
!    this%order = order
!  end subroutine node_array_create

!  subroutine node_array_fill ( this )
!    implicit none
!    class(node_array_t), intent(inout) :: this
!    type(node_iterator_t) :: node_iterator 
!    integer(ip) :: i, c, m_ijk, current
!    current = ISHFT(ISHFT(1_ip,DIM)-1,DIM)
!    node_iterator = this%create_node_iterator( current, own_boundary = .true. )
!    call node_iterator%init()
!    c = 0
!    m_ijk = 0
!    do while (.not. node_iterator%has_finished() )
!       c = c+1
!       m_ijk = max( m_ijk, node_iterator%current() )
!       call node_iterator%next()       
!    end do
!    this%number_nodes = c
!    call memalloc ( c, this%node_array, __FILE__, __LINE__ )
!    call memalloc ( m_ijk+1, this%ijk_to_index, __FILE__, __LINE__, lb1 = 0 )
!    call node_iterator%init()
!    c = 0
!    do while (.not. node_iterator%has_finished() )
!       c = c+1
!       this%node_array(c) =  node_iterator%current()
!       call node_iterator%next()         
!    end do
!    ! Sort nodes based on its bit-based index
!    call sort( c, this%node_array )    
!    write(*,*) 'node_array',this%node_array
!    do i =1,c
!       write(*,*) 'object(',i,') :'
!       write(*,'(B32)') this%node_array(i)
!       this%ijk_to_index( this%node_array(i) ) = i
!    end do
!  end subroutine node_array_fill

!  subroutine node_array_print( this )
!    implicit none
!    class(node_array_t), intent(inout) :: this
!    integer(ip) :: i
!    write(*,*) '************NODE ARRAY***********'
!    write(*,*) 'number of nodes = ', this%number_nodes
!    write(*,*) 'object_tree_topology = '
!    write(*,'(B32)') this%object_tree%topology
!    write(*,*) 'nodes'
!    do i = 1, this%number_nodes
!       write(*,*) 'NODE: ',this%node_array(i)
!       write(*,'(B32)') this%node_array(i)
!    end do
!    write(*,*) 'ijk_to_index: ',this%ijk_to_index
!  end subroutine node_array_print

!end module abstract_reference_element_names


program abstract_reference_element
  use serial_names
  !use abstract_reference_element_names
  implicit none
  !type(geometry_tree_t) :: object_tree
  !type(children_iterator_t) :: children_iterator
  !type(node_array_t) :: node_array
  !type(node_iterator_t) :: node_iterator
  !integer(ip), parameter :: DIM = 2
  !integer(ip) :: root, children, node, c, i

  call fempar_init()
  
  ! Assuming we are in a reference_fe, and we want to fill
  
  ! interior_nodes_vef
  
  
  
  
  
  
  
  
  
  
  
  
  
  !call object_tree%create( topology )  
  !call node_array%create( object_tree, order )
  !call node_array%fill()
  !!
  !ref_fe%number_nodes = node_array%number_nodes
  !ref_fe%number_vefs = object_tree%number_vefs
  !!  
  !ref_fe%number_vef_dimensions = 0
  !do i = 1,object_tree%number_objects
  !   obj_dim = object_dimension( object_tree%object_array(i) )
  !   ref_fe%number_vefs_dimension( obj_dim+1 ) = ref_fe%number_vefs_dimension( obj_dim+1 ) + 1
  !end do  
  !! interior nodes vef
  !do i = 1,object_tree%number_objects
  !   node_iterator = node_array%create_node_iterator( object_tree%object_array(i), own_boundary = .false. )
  !   do while (.not. node_iterator%has_finished() )
  !      node = node_iterator%current()
  !      ! XXX
  !      call node_iterator%next()
  !   end do
  !   ! XXX
  !end do
  !! nodes_vef
  !do i = 1,object_tree%number_objects
  !   node_iterator = node_array%create_node_iterator( object_tree%object_array(i), own_boundary = .true. )
  !   do while (.not. node_iterator%has_finished() )
  !      node = node_iterator%current()
  !      ! XXX
  !      call node_iterator%next()
  !   end do
  !end do
  !! vertices_vef
  !call vertex_array%create( object_tree, order = 1 )
  !call vertex_array%fill()
  !do i = 1,object_tree%number_objects
  !   node_iterator = vertex_array%create_node_iterator( object_tree%object_array(i), own_boundary = .true. )
  !   do while (.not. node_iterator%has_finished() )
  !      node = node_iterator%current()
  !      ! XXX
  !      call node_iterator%next()
  !   end do
  !end do
  !! vefs_vef (NEW DEFINITION)
  !do i = 1,object_tree%number_objects
  !   children_iterator = this%create_children_iterator( object_tree%object_array(i) )
  !   do while (.not. children_iterator%has_finished() )
  !      children = children_iterator%current()
  !      ! XXX
  !      call children_iterator%next()
  !   end do
  !end do
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ! Easy part...
  ! own_nodes_vef
  ! face_integration_coupling_nodes_face
  
  
  !!call object_tree%print()
  !call node_array%print()

  !do i=1,object_tree%number_objects
  !   write(*,*) '***************** OBJECT :',i, ' IJK: *****************'
  !   write(*,'(B32)') object_tree%object_array(i)
  !   node_iterator = node_array%create_node_iterator( object_tree%object_array(i), own_boundary = .true. )
  !   call node_iterator%init()  
  !   write(*,*) 'OBJECT NODES:'
  !   do while (.not. node_iterator%has_finished() )
  !      node = node_iterator%current()
  !      write(*,*) 'NODE COORDINATES: ', node_iterator%coordinate
  !      write(*,*) 'NODE: ', node_array%ijk_to_index(node)
  !      write(*,*) 'NODE IJK: ' 
  !      write(*,'(B32)') node
  !      
  !      call node_iterator%next()
  !      !   node = node_iterator%current()
  !      !   c = c+1
  !   end do
  !end do






  !children = 56
  !children = 49
  !children = 12
  !children_iterator = object_tree%create_children_iterator(children)


  ! Node numbering






  ! Loop over children
  !call children_iterator%init()
  !do while (.not. children_iterator%has_finished() )
  !   children = children_iterator%current()
  !   write(*,*) 'current children: '
  !   write(*,'(B32)') children
  ! Put children of child (recursive call)
  !call this%fill_tree( children )




  !node_iterator = object_tree%create_node_iterator( children, order = 3, own_boundary = .false. )
  !! Loop over children
  !call node_iterator%init()
  !c = 0
  !do while (.not. node_iterator%has_finished() )
  !   node = node_iterator%current()
  !   write(*,*) 'current displacement: ', node_iterator%displacement
  !   write(*,*) 'current node: ', node
  !   call node_iterator%next()
  !   node = node_iterator%current()
  !   c = c+1
  !end do
  !write(*,*) 'Number of nodes:',c

  ! Move to the next child
  !   call children_iterator%next()
  !end do




  call fempar_finalize()





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



  !write(*,*) 'DIMENSION'
  !write(*,'(B32)') DIM
end program abstract_reference_element
