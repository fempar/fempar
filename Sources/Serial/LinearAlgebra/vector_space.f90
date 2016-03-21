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
module vector_space_names
  use types_names
  use vector_names
  implicit none

  private
  
  ! type(vector_space_t) is a concrete data type which represents a vector 
  ! space in its strict mathematical sense. This data type follows two different 
  ! creational patterns:
  ! 1. PROTOTYPE. It aggregates a polymorphic class(vector_t) instance,
  !    that acts as a prototypical instance.
  ! 2. FACTORY_METHOD. It provides a TBP, called create_vector, which
  !    creates a class(vector_t) polymorphic instance by cloning the
  !    prototypical instance. create_vector is the factory method.
  type:: vector_space_t
    private
    class(vector_t), allocatable :: vector
   contains
     procedure :: create        => vector_space_create
     procedure :: is_created    => vector_space_is_created
     procedure :: create_vector => vector_space_create_vector
     procedure :: belongs_to    => vector_space_belongs_to
     procedure :: equal_to      => vector_space_equal_to
     procedure :: clone         => vector_space_clone 
     procedure :: free          => vector_space_free
					procedure :: get_number_blocks
  end type vector_space_t

  public :: vector_space_t

contains  

     subroutine vector_space_create(this,vector)
       implicit none
       class(vector_space_t), intent(inout) :: this
       class(vector_t)      , intent(in)    :: vector
       call this%free()
       allocate(this%vector, mold=vector)
       call this%vector%clone(vector)
     end subroutine vector_space_create


     function vector_space_is_created(this) result(is_created)
       implicit none
       class(vector_space_t), intent(in) :: this
       logical                           :: is_created
       is_created = allocated(this%vector)
     end function vector_space_is_created

     ! This method selects the dynamic type of class(vector_t) and allocates it, 
     ! following the FACTORY METHOD OOD pattern. As the method is responsible for
     ! selecting the dynamic type of vector, the polymorphic entity should have the
     ! pointer or allocatable attribute. I decided to have the allocatable attribute instead of 
     ! pointer to be aligned with the way that abstract data calculus among class(vector_t) polymorphic 
     ! instances is managed. I do not see how this method can be a function instead of a subroutine. 
     ! If this would be a function, then the polymorphic vector instance returned by the function would be temporal, so that
     ! the client should find a way to transfer the temporary returned allocatable to 
     ! the permanent one handled by the client. This cannot be achieved with the move_alloc
     ! intrinsic (https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/495399). 
     ! While the allocate(permanent,source=f()) construct is legal in this context, this would require an additional copy. 
     ! The copy should be performed with the overloaded assignment, instead of the implicit one. I do not know whether 
     ! the allocate(*,source=*) calls overloaded assignment or not. 
     subroutine vector_space_create_vector(this,vector) 
       implicit none
       class(vector_space_t)             , intent(in)    :: this
       class(vector_t)      , allocatable, intent(inout) :: vector
       
       if (allocated(vector)) then
          call vector%free()
          deallocate(vector)
       end if
       
       allocate(vector, mold=this%vector)
       call vector%clone(this%vector)
       call vector%allocate()
     end subroutine vector_space_create_vector
	
     ! Determines the dynamic type of the pointer result to match that of this 
     ! (when allocating the pointer result), copying the contents of this into the result 
     ! of the function. This deferred TBP only has sense in case of an aggregation/composition
     ! relantionship among operator_t and vector_space_t
     subroutine vector_space_clone(this,vector_space)
       implicit none
       class(vector_space_t), intent(in)    :: this
       type(vector_space_t), intent(inout)  :: vector_space
       call vector_space%free()
       allocate (vector_space%vector, mold=this%vector)
       call vector_space%vector%clone(this%vector)
     end subroutine vector_space_clone

     ! Determine whether the input vector_space is equal to the passed-object dummy argument
     ! vector space this.
     function vector_space_equal_to(this, vector_space)
       implicit none
       class(vector_space_t), intent(in) :: this
       class(vector_space_t), intent(in) :: vector_space
       logical :: vector_space_equal_to
       vector_space_equal_to = this%vector%same_vector_space(vector_space%vector)
     end function vector_space_equal_to

     ! Returns true if vector belongs to this and false if it does not
     function vector_space_belongs_to(this,vector)
       implicit none
       class(vector_space_t), intent(in) :: this
       class(vector_t)      , intent(in) :: vector
       logical vector_space_belongs_to
       vector_space_belongs_to = this%vector%same_vector_space(vector)
     end function vector_space_belongs_to
     
     subroutine vector_space_free(this)
       implicit none
       class(vector_space_t), intent(inout) :: this
       if (allocated(this%vector)) then
         call this%vector%free()
         deallocate(this%vector)
       end if  
     end subroutine vector_space_free
					
			function get_number_blocks(this) result(res)
     implicit none 
     class(vector_space_t), intent(in) :: this
     integer(ip) :: res
			  res = this%vector%get_number_blocks()
   end function get_number_blocks
     
end module vector_space_names
