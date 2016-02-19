module list_types_names

USE types_names, ONLY: ip
USE memor_names

implicit none

    type list_t
    !-----------------------------------------------------------------
    !< List_t derived type
    !----------------------------------------------------------------- 
        integer(ip)              :: n
        integer(ip), allocatable :: p(:) 
        integer(ip), allocatable :: l(:) 
    contains
        procedure, private :: list_assign
        procedure          :: create                     => list_create
        procedure          :: allocate_list_from_pointer => list_allocate_list_from_p
        procedure          :: free                       => list_free
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
		

contains

!--------------------------------------------------------------------- -----------------------------------------------------------
!< List_t derived type TBP's
!--------------------------------------------------------------------- -----------------------------------------------------------

    subroutine list_free( this )
        class(list_t), intent(inout) :: this
        if(allocated(this%p)) call memfree(this%p, __FILE__, __LINE__)
        if(allocated(this%l)) call memfree(this%l, __FILE__, __LINE__)
        this%n = 0
    end subroutine list_free
  
    subroutine list_create( this, n )
        class(list_t), intent(inout) :: this
        integer(ip)  ,   intent(in)  :: n
        call this%free()
        this%n = n
        call memalloc(this%n+1, this%p, __FILE__, __LINE__)
        this%p = 0
    end subroutine list_create  
  
    subroutine list_allocate_list_from_p( this )
        class(list_t), intent(inout) :: this
        call memalloc(this%p(this%n+1)-1, this%l, __FILE__, __LINE__)
        this%l = 0
    end subroutine list_allocate_list_from_p  
  
    subroutine list_assign(this,list) 
        class(list_t), intent(inout):: this
        type(list_t), intent(in):: list
        call this%free()
        call this%create(n=list%n)
        this%p(:) = list%p(:)
        call this%allocate_list_from_pointer()
        this%l(:) = list%l(:)
    end subroutine list_assign
  
    subroutine list_print( this, lunou )
        class(list_t)    , intent(in) :: this
        integer(ip)      , intent(in) :: lunou
        integer(ip) :: i
        write(lunou,*) '****PRINT LIST****'
        write(lunou,*) 'size total list:',this%n
        do i = 1,this%n
            write(lunou,*) 'l(',i,')',this%l(this%p(i):this%p(i+1)-1)
        end do
        write(lunou,*) '****END PRINT LIST****'
    end subroutine list_print

!--------------------------------------------------------------------- -----------------------------------------------------------
!< List_t derived type TBP's
!--------------------------------------------------------------------- -----------------------------------------------------------

    subroutine list_2d_free( this )
        class(list_2d_t), intent(inout) :: this
        if(allocated(this%p)) call memfree(this%p, __FILE__, __LINE__)
        if(allocated(this%l)) call memfree(this%l, __FILE__, __LINE__)
        this%n1 = 0
        this%n2 = 0
    end subroutine list_2d_free  

    subroutine list_2d_create( this, n1, n2 )
        class(list_2d_t), intent(inout) :: this
        integer(ip),      intent(in)    :: n1
        integer(ip),      intent(in)    :: n2
	    call this%free()
        this%n1 = n1
        this%n2 = n2
        call memalloc(this%n1+1, this%p, __FILE__, __LINE__)
        this%p = 0
    end subroutine list_2d_create  

    subroutine list_2d_allocate_list_from_p( this )
        class(list_2d_t), intent(inout) :: this
        call memalloc(this%p(this%n1+1)-1, this%n2, this%l, __FILE__, __LINE__)
        this%l = 0
    end subroutine list_2d_allocate_list_from_p  
  
    subroutine list_2d_assign(this,list) 
        class(list_2d_t), intent(inout):: this
        type(list_2d_t), intent(in):: list
        call this%free()
        call this%create(n1=list%n1, n2=list%n2)
        this%p(:) = list%p(:)
        call this%allocate_list_from_pointer()
        this%l(:,:) = list%l(:,:)
    end subroutine list_2d_assign  
  
    subroutine list_2d_print( this, lunou )
        class(list_2d_t), intent(in) :: this
        integer(ip) ,     intent(in) :: lunou
        integer(ip) :: i,j
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

end module list_types_names
