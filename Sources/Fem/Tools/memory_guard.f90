module memory_guard_names
  use types
  implicit none
  private
  
  ! Object-Oriented design pattern taken from the following paper:
  ! D. W. I. Rouson, J. Xia, X. Xu
  ! Object construction and destruction design patterns in Fortran 2003
  ! Procedia Computer Science 1 (2012), 1495--1504
  type, abstract :: memory_guard
     private
     integer(ip), pointer :: temporary =>null() ! Null marks non-temporary data
   contains
     procedure  :: SetTemp   ! Mark object as temporary
     procedure  :: GuardTemp ! Increment the depth count
     procedure  :: CleanTemp ! Decrement depth count/Free memory if 1
     procedure(free_interface), deferred :: free
  end type memory_guard
  
  abstract interface
     subroutine free_interface ( this )
       import :: memory_guard
       class(memory_guard), intent(inout) :: this
     end subroutine free_interface
  end interface

  public :: memory_guard

contains

  subroutine SetTemp ( this )
    implicit none
    class(memory_guard), intent(inout) :: this
    if ( .not. associated(this%temporary)) allocate(this%temporary)
    this%temporary = 1 
  end subroutine SetTemp

  subroutine GuardTemp ( this )
    implicit none
    class(memory_guard), intent(in) :: this
    if (associated(this%temporary)) this%temporary = this%temporary + 1 
  end subroutine GuardTemp

  subroutine CleanTemp ( this )
    implicit none
    class(memory_guard) :: this
    if (associated(this%temporary)) then
       if (this%temporary > 1) this%temporary = this%temporary - 1
       if (this%temporary == 1) then
          call this%free()
          deallocate(this%temporary)
       end if
    end if
  end subroutine CleanTemp
  
end module memory_guard_names
