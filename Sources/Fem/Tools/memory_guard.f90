module memory_guard_names
use types_names
  implicit none
  private
  
  ! Object-Oriented design pattern taken from the following paper:
  ! D. W. I. Rouson, J. Xia, X. Xu
  ! Object construction and destruction design patterns in Fortran 2003
  ! Procedia Computer Science 1 (2012), 1495--1504
  type, abstract :: memory_guard_t
     private
     integer(ip), pointer :: temporary =>null() ! Null marks non-temporary data
   contains
     procedure  :: SetTemp   ! Mark object as temporary
     procedure  :: GuardTemp ! Increment the depth count
     procedure  :: CleanTemp ! Decrement depth count/Free memory if 1
     procedure  :: IsTemp
     procedure  :: GetTemp
     procedure  :: PrintTemp ! Print value
     procedure(free_interface), deferred :: free
  end type memory_guard_t
  
  abstract interface
     subroutine free_interface ( this )
       import :: memory_guard_t
       class(memory_guard_t), intent(inout) :: this
     end subroutine free_interface
  end interface

  public :: memory_guard_t

contains

  subroutine SetTemp ( this )
    implicit none
    class(memory_guard_t), intent(inout) :: this
    if ( .not. associated(this%temporary)) allocate(this%temporary)
    this%temporary = 1 
  end subroutine SetTemp

  subroutine GuardTemp ( this )
    implicit none
    class(memory_guard_t), intent(in) :: this
    if (associated(this%temporary)) this%temporary = this%temporary + 1 
  end subroutine GuardTemp

  subroutine CleanTemp ( this )
    implicit none
    class(memory_guard_t) :: this
    if (associated(this%temporary)) then
       if (this%temporary > 1) this%temporary = this%temporary - 1
       if (this%temporary == 1) then
          call this%free()
          deallocate(this%temporary)
       end if
    end if
  end subroutine CleanTemp

  function IsTemp( this ) 
    implicit none
   class(memory_guard_t) :: this
   logical :: IsTemp
   IsTemp = associated(this%temporary)
 end function IsTemp
 
 function GetTemp( this ) 
   implicit none
   class(memory_guard_t) :: this
   integer(ip) :: GetTemp
   if( associated(this%temporary)) then
      GetTemp=this%temporary
   else
      GetTemp=-1
   end if
 end function GetTemp
 
 subroutine PrintTemp ( this )
   implicit none
   class(memory_guard_t), intent(in) :: this
   if ( associated(this%temporary)) then
      write(*,*) 'It is temporary with flag', this%temporary
   else
      write(*,*) 'It is permanent'
   end if
 end subroutine PrintTemp
 
end module memory_guard_names
