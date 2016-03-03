module list_types_names

USE types_names, ONLY: ip
USE memor_names

implicit none

# include "debug.i90"

    !-----------------------------------------------------------------
    !< Transition diagram STATES 
    !----------------------------------------------------------------- 
    integer(ip), parameter :: LIST_STATE_START          = 0
    integer(ip), parameter :: LIST_STATE_CREATED        = 1
    integer(ip), parameter :: LIST_STATE_HEADER_BUILT   = 2
    integer(ip), parameter :: LIST_STATE_LIST_ALLOCATED = 3
    !-----------------------------------------------------------------
    ! Input State         | Action                | Output State 
    !-----------------------------------------------------------------
    ! START               | free()                | START
    ! START               | create()              | CREATED
    !
    ! CREATED             | free()                | START
    ! CREATED             | create()              | CREATED
    ! CREATED             | calculate_header()    | HEADER_BUILT
    !
    ! HEADER_BUILT        | free()                | START
    ! HEADER_BUILT        | create()              | CREATED
    ! HEADER_BUILT        | allocate_list_from_p()| LIST_ALLOCATED
    !
    ! LIST_ALLOCATED      | free()                | START
    ! LIST_ALLOCATED      | create()              | CREATED
    !-----------------------------------------------------------------

    type list_t
    !-----------------------------------------------------------------
    !< List_t derived type
    !----------------------------------------------------------------- 
        integer(ip), private     :: state = LIST_STATE_START
        integer(ip)              :: n
        integer(ip), allocatable :: p(:) 
        integer(ip), allocatable :: l(:) 
    contains
        procedure, private :: list_assign
        procedure, private :: list_get_full_list_iterator
        procedure, private :: list_get_list_range_iterator
        procedure, private :: list_get_list_index_iterator
        procedure          :: create                     => list_create
        procedure          :: get_num_pointers           => list_get_num_pointers
        procedure          :: get_size                   => list_get_size
        procedure          :: sum_to_pointer_index       => list_sum_to_pointer_index
        procedure          :: calculate_header           => list_calculate_header
        procedure          :: allocate_list_from_pointer => list_allocate_list_from_p
        procedure          :: free                       => list_free
        generic            :: get_iterator               => list_get_full_list_iterator,  &
                                                            list_get_list_range_iterator, &
                                                            list_get_list_index_iterator
        generic            :: assignment(=)              => list_assign
        procedure          :: print                      => list_print
    end type list_t

    type list_2d_t
    !-----------------------------------------------------------------
    !< List_2d_t
    !----------------------------------------------------------------- 
        integer(ip)              :: n1
        integer(ip)              :: n2
        integer(ip), allocatable :: p(:) 
        integer(ip), allocatable :: l(:,:) 
    contains
        procedure, private :: list_2d_assign
        procedure          :: create                     => list_2d_create
        procedure          :: allocate_list_from_pointer => list_2d_allocate_list_from_p
        procedure          :: free                       => list_2d_free
        generic            :: assignment(=)              => list_2d_assign
        procedure          :: print                      => list_2d_print
    end type list_2d_t

    type list_pointer_t
    !-----------------------------------------------------------------
    !< List_pointer_t derived type
    !----------------------------------------------------------------- 
        type(list_t), pointer :: p => NULL()
    end type list_pointer_t


    type list_iterator_t
    private
    !-----------------------------------------------------------------
    !< List_t derived type
    !----------------------------------------------------------------- 
        type(list_pointer_t)  :: list 
        integer(ip)           :: first
        integer(ip)           :: current
        integer(ip)           :: last
    contains
    private
        procedure         :: init                        => list_iterator_init
        procedure         :: is_in_range                 => list_iterator_is_in_range
        procedure, public :: is_lower_bound              => list_iterator_is_lower_bound
        procedure, public :: is_upper_bound              => list_iterator_is_upper_bound
        procedure, public :: has_previous                => list_iterator_has_previous
        procedure, public :: has_next                    => list_iterator_has_next
        procedure, public :: get_size                    => list_iterator_get_size
        procedure, public :: get_distance_to_lower_bound => list_iterator_get_distance_to_lower_bound
        procedure, public :: get_distance_to_upper_bound => list_iterator_get_distance_to_upper_bound
        procedure, public :: begin                       => list_iterator_begin
        procedure, public :: end                         => list_iterator_end
        procedure, public :: previous                    => list_iterator_previous
        procedure, public :: next                        => list_iterator_next
        procedure, public :: get_current                 => list_iterator_get_current
        procedure, public :: set_current                 => list_iterator_set_current
        procedure, public :: reach_from_current          => list_iterator_reach_from_current
    end type list_iterator_t

contains

!--------------------------------------------------------------------- -----------------------------------------------------------
!< List_t derived type TBP's
!--------------------------------------------------------------------- -----------------------------------------------------------

    subroutine list_free( this )
    !-----------------------------------------------------------------
    !< Free list_T derived type
    !----------------------------------------------------------------- 
        class(list_t), intent(inout) :: this
    !----------------------------------------------------------------- 
        if(allocated(this%p)) call memfree(this%p, __FILE__, __LINE__)
        if(allocated(this%l)) call memfree(this%l, __FILE__, __LINE__)
        this%n = 0
        this%n = LIST_STATE_START
    end subroutine list_free


    subroutine list_create( this, n )
    !-----------------------------------------------------------------
    !< Allocate pointer and initialize it to zero
    !----------------------------------------------------------------- 
        class(list_t), intent(inout) :: this
        integer(ip),   intent(in)    :: n
    !----------------------------------------------------------------- 
        call this%free()
        this%n = n
        call memalloc(this%n+1, this%p, __FILE__, __LINE__)
        this%p = 0
        this%state = LIST_STATE_CREATED
    end subroutine list_create  


    subroutine list_sum_to_pointer_index( this, index, value )
    !-----------------------------------------------------------------
    !< Sum a value given a pointer index
    !----------------------------------------------------------------- 
        class(list_t), intent(inout) :: this
        integer(ip),   intent(in)    :: index
        integer(ip),   intent(in)    :: value
    !----------------------------------------------------------------- 
!        assert(this%state == LIST_STATE_CREATED)
        this%p(index+1) = this%p(index+1) + value
    end subroutine list_sum_to_pointer_index


    subroutine list_calculate_header( this)
    !-----------------------------------------------------------------
    !< Transforms the pointer from counter to header
    !----------------------------------------------------------------- 
        class(list_t), intent(inout) :: this
        integer(ip)                  :: i
    !----------------------------------------------------------------- 
!        assert(this%state == LIST_STATE_CREATED)
        this%p(1) = 1
        do i=1, this%n
            this%p(i+1) = this%p(i+1) + this%p(i)
        enddo
        this%state = LIST_STATE_HEADER_BUILT
    end subroutine list_calculate_header

  
    subroutine list_allocate_list_from_p( this )
    !----------------------------------------------------------------- 
    !< Allocate list_t%l using the pointer info
    !----------------------------------------------------------------- 
        class(list_t), intent(inout) :: this
    !----------------------------------------------------------------- 
!        assert(this%state == LIST_STATE_HEADER_BUILT)
        call memalloc(this%p(this%n+1)-1, this%l, __FILE__, __LINE__)
        this%l = 0
        this%state = LIST_STATE_LIST_ALLOCATED
    end subroutine list_allocate_list_from_p  


    function list_get_full_list_iterator( this ) result(list_iterator)
    !----------------------------------------------------------------- 
    !< Allocate list_t%l using the pointer info
    !----------------------------------------------------------------- 
        class(list_t), target, intent(inout) :: this
        type(list_iterator_t)                :: list_iterator
    !----------------------------------------------------------------- 
!        assert(this%state == LIST_STATE_LIST_STATE_LIST_ALLOCATED)
        call list_iterator%init(list=this, first=this%p(1), current=this%p(1), last=this%p(this%n+1)-1)
    end function list_get_full_list_iterator


    function list_get_list_range_iterator( this, start, end ) result(list_iterator)
    !----------------------------------------------------------------- 
    !< Allocate list_t%l using the pointer info
    !----------------------------------------------------------------- 
        class(list_t), target, intent(inout) :: this
        integer(ip),           intent(in)    :: start
        integer(ip),           intent(in)    :: end
        type(list_iterator_t)                :: list_iterator
    !----------------------------------------------------------------- 
!        assert(this%state == LIST_STATE_LIST_STATE_LIST_ALLOCATED)
        assert( (start <= end) .and. (start>=1) .and. (end<= this%n) )
        call list_iterator%init(list=this, first=this%p(start), current=this%p(start), last=this%p(end+1)-1)
    end function list_get_list_range_iterator


    function list_get_list_index_iterator( this, index ) result(list_iterator)
    !----------------------------------------------------------------- 
    !< Allocate list_t%l using the pointer info
    !----------------------------------------------------------------- 
        class(list_t), target, intent(inout) :: this
        integer(ip),           intent(in)    :: index
        type(list_iterator_t)                :: list_iterator
    !----------------------------------------------------------------- 
!        assert(this%state == LIST_STATE_LIST_STATE_LIST_ALLOCATED)
        assert( (index>=1) .and. (index<= this%n) )
        call list_iterator%init(list=this, first=this%p(index), current=this%p(index), last=this%p(index+1)-1)
    end function list_get_list_index_iterator


    subroutine list_assign(this,list) 
    !----------------------------------------------------------------- 
    !< Assign operator between list_t derived types
    !----------------------------------------------------------------- 
        class(list_t), intent(inout):: this
        type(list_t), intent(in):: list
    !----------------------------------------------------------------- 
!        assert(list%state == LIST_STATE_LIST_ALLOCATED)
        call this%free()
        call this%create(n=list%n)
        this%p(:) = list%p(:)
        this%state = LIST_STATE_HEADER_BUILT
        call this%allocate_list_from_pointer()
        this%l(:) = list%l(:)
    end subroutine list_assign


    function list_get_num_pointers( this ) result(num_pointers)
    !-----------------------------------------------------------------
    !< Return the number of pointers
    !----------------------------------------------------------------- 
        class(list_t), intent(inout) :: this
        integer(ip)                  :: num_pointers
    !----------------------------------------------------------------- 
!        assert(this%state == LIST_STATE_CREATED)
        num_pointers = this%n
    end function list_get_num_pointers


    function list_get_size( this ) result(list_size)
    !-----------------------------------------------------------------
    !< Return the size of the list
    !----------------------------------------------------------------- 
        class(list_t), intent(inout) :: this
        integer(ip)                  :: list_size
    !----------------------------------------------------------------- 
!        assert(this%state == LIST_STATE_LIST_ALLOCATED)
        list_size = size(this%l)
    end function list_get_size


    subroutine list_print( this, lunou )
    !----------------------------------------------------------------- 
    !< Print the list content
    !----------------------------------------------------------------- 
        class(list_t)    , intent(in) :: this
        integer(ip)      , intent(in) :: lunou
        integer(ip) :: i
    !----------------------------------------------------------------- 
        write(lunou,*) '****PRINT LIST****'
        write(lunou,*) 'size total list:',this%n
        do i = 1,this%n
            write(lunou,*) 'l(',i,')',this%l(this%p(i):this%p(i+1)-1)
        end do
        write(lunou,*) '****END PRINT LIST****'
    end subroutine list_print

!--------------------------------------------------------------------- -----------------------------------------------------------
!< List_2d_t derived type TBP's
!--------------------------------------------------------------------- -----------------------------------------------------------

    subroutine list_2d_free( this )
    !----------------------------------------------------------------- 
    !< Free the list_2d_t derived type
    !----------------------------------------------------------------- 
        class(list_2d_t), intent(inout) :: this
    !----------------------------------------------------------------- 
        if(allocated(this%p)) call memfree(this%p, __FILE__, __LINE__)
        if(allocated(this%l)) call memfree(this%l, __FILE__, __LINE__)
        this%n1 = 0
        this%n2 = 0
    end subroutine list_2d_free  


    subroutine list_2d_create( this, n1, n2 )
    !----------------------------------------------------------------- 
    !< Allocate pointer and initialize it to zero and set the second
    !< dimension of %l
    !-----------------------------------------------------------------
        class(list_2d_t), intent(inout) :: this
        integer(ip),      intent(in)    :: n1
        integer(ip),      intent(in)    :: n2
    !-----------------------------------------------------------------
	    call this%free()
        this%n1 = n1
        this%n2 = n2
        call memalloc(this%n1+1, this%p, __FILE__, __LINE__)
        this%p = 0
    end subroutine list_2d_create  


    subroutine list_2d_allocate_list_from_p( this )
    !----------------------------------------------------------------- 
    !< Allocate the list from pointer info
    !----------------------------------------------------------------- 
        class(list_2d_t), intent(inout) :: this
    !----------------------------------------------------------------- 
        call memalloc(this%p(this%n1+1)-1, this%n2, this%l, __FILE__, __LINE__)
        this%l = 0
    end subroutine list_2d_allocate_list_from_p  


    subroutine list_2d_assign(this,list) 
    !----------------------------------------------------------------- 
    !< Assign operator between list_2d_t derived types
    !----------------------------------------------------------------- 
        class(list_2d_t), intent(inout):: this
        type(list_2d_t), intent(in):: list
    !----------------------------------------------------------------- 
        call this%free()
        call this%create(n1=list%n1, n2=list%n2)
        this%p(:) = list%p(:)
        call this%allocate_list_from_pointer()
        this%l(:,:) = list%l(:,:)
    end subroutine list_2d_assign  

  
    subroutine list_2d_print( this, lunou )
    !----------------------------------------------------------------- 
    !< Print the list content
    !----------------------------------------------------------------- 
        class(list_2d_t), intent(in) :: this
        integer(ip) ,     intent(in) :: lunou
        integer(ip) :: i,j
    !----------------------------------------------------------------- 
        write(lunou,*) '****PRINT LIST_2D****'
        write(lunou,*) 'size total list:',this%n1
        write(lunou,*) 'size components:',this%n2
        do j = 1,this%n2
            do i = 1,this%n1
                write(lunou,*) 'l(',i,',',j,')',this%l(this%p(i):this%p(i+1)-1,j)
            end do
        end do
        write(lunou,*) '****END PRINT LIST_2D****'
    end subroutine list_2d_print

!--------------------------------------------------------------------- -----------------------------------------------------------
!< List_iterator_t derived type TBP's
!--------------------------------------------------------------------- -----------------------------------------------------------

    subroutine list_iterator_init(this, list, first, current, last)
    !----------------------------------------------------------------- 
    !< A list_T initializes an iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        type(list_t), target,   intent(in)    :: list
        integer(ip),            intent(in)    :: first
        integer(ip),            intent(in)    :: current
        integer(ip),            intent(in)    :: last
    !----------------------------------------------------------------- 
        this%list%p  => list
        this%first   = first
        this%current = current
        this%last    = last
    end subroutine list_iterator_init


    function list_iterator_is_in_range(this, index) result(in_range)
    !----------------------------------------------------------------- 
    !< Check if a given index is in iterator range
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        integer(ip),            intent(in)    :: index
        logical                               :: in_range
    !----------------------------------------------------------------- 
        in_range = ( (index>=this%first) .and. (index<=this%last) )
    end function list_iterator_is_in_range


    function list_iterator_get_size(this) result(size)
    !----------------------------------------------------------------- 
    !< Returns the total size of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        integer(ip)                           :: size
    !----------------------------------------------------------------- 
        size = this%last-this%first+1
    end function list_iterator_get_size


    function list_iterator_get_distance_to_lower_bound(this) result(distance)
    !----------------------------------------------------------------- 
    !< Returns the distance from the current position (included) to the
    !< lower bound of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        integer(ip)                           :: distance
    !----------------------------------------------------------------- 
        distance = this%current-this%first+1
    end function list_iterator_get_distance_to_lower_bound


    function list_iterator_get_distance_to_upper_bound(this) result(distance)
    !----------------------------------------------------------------- 
    !< Returns the distance from the current position (included) to the
    !< upper bound of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        integer(ip)                           :: distance
    !----------------------------------------------------------------- 
        distance = this%last-this%current+1
    end function list_iterator_get_distance_to_upper_bound


    subroutine list_iterator_begin(this)
    !----------------------------------------------------------------- 
    !< Jump to first position of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
    !----------------------------------------------------------------- 
        this%current = this%first
    end subroutine list_iterator_begin


    subroutine list_iterator_end(this)
    !----------------------------------------------------------------- 
    !< Jump to last position of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
    !----------------------------------------------------------------- 
        this%current = this%last
    end subroutine list_iterator_end


    function list_iterator_is_lower_bound(this) result(is_lower_bound)
    !----------------------------------------------------------------- 
    !< Check if the current component is a lower bound of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        logical                               :: is_lower_bound
    !----------------------------------------------------------------- 
        is_lower_bound = (this%current < this%first)
    end function list_iterator_is_lower_bound


    function list_iterator_is_upper_bound(this) result(is_upper_bound)
    !----------------------------------------------------------------- 
    !< Check if the current component is a upper bound of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        logical                               :: is_upper_bound
    !----------------------------------------------------------------- 
        is_upper_bound = (this%current > this%last)
    end function list_iterator_is_upper_bound


    function list_iterator_has_next(this) result(has_next)
    !----------------------------------------------------------------- 
    !< Check if the current component has next
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        logical                               :: has_next
    !----------------------------------------------------------------- 
        has_next = (this%current < this%last)
    end function list_iterator_has_next


    function list_iterator_has_previous(this) result(has_previous)
    !----------------------------------------------------------------- 
    !< Check if the current component has next
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        logical                               :: has_previous
    !----------------------------------------------------------------- 
        has_previous = (this%current > this%first)
    end function list_iterator_has_previous


    function list_iterator_get_current(this) result(component)
    !----------------------------------------------------------------- 
    !< Returns the current component of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        integer(ip)                           :: component
    !----------------------------------------------------------------- 
        component = this%list%p%l(this%current)
    end function list_iterator_get_current


    subroutine list_iterator_set_current(this, value)
    !----------------------------------------------------------------- 
    !< Returns the current component of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        integer(ip),            intent(in)    :: value
    !----------------------------------------------------------------- 
        this%list%p%l(this%current) = value
    end subroutine list_iterator_set_current


    function list_iterator_reach_from_current(this, offset) result(component)
    !----------------------------------------------------------------- 
    !< Returns the current + offset component of the iterator
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
        integer(ip),            intent(in)    :: offset
        integer(ip)                           :: component
    !----------------------------------------------------------------- 
        assert( this%is_in_range(this%current+offset) )
        component = this%list%p%l(this%current+offset)
    end function list_iterator_reach_from_current


    subroutine list_iterator_previous(this)
    !----------------------------------------------------------------- 
    !< Move the cursor of the iterator to the previous component
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
    !----------------------------------------------------------------- 
        if(.not. this%is_lower_bound()) this%current = this%current-1
    end subroutine list_iterator_previous


    subroutine list_iterator_next(this)
    !----------------------------------------------------------------- 
    !< Move the cursor of the iterator to the next component
    !----------------------------------------------------------------- 
        class(list_iterator_t), intent(inout) :: this
    !----------------------------------------------------------------- 
        if(.not. this%is_upper_bound()) this%current = this%current+1
    end subroutine list_iterator_next


end module list_types_names
