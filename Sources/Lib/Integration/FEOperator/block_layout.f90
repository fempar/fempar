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
module block_layout_names
  use types_names
  use memor_names

  implicit none
# include "debug.i90"
  private

  type :: block_layout_t
   private
   integer(ip)               :: num_blocks = -1
   integer(ip)               :: num_fields = -1
   integer(ip), allocatable  :: field_id_to_block_id(:)
   integer(ip), allocatable  :: num_dofs_x_block(:)
   logical    , allocatable  :: field_coupling(:,:)
   logical    , allocatable  :: block_coupling(:,:)
  contains
     procedure, non_overridable :: create                   => block_layout_create
     procedure, non_overridable :: free                     => block_layout_free
     procedure, non_overridable :: clear_num_dofs_x_block   => block_layout_clear_num_dofs_x_block
     procedure, non_overridable :: get_num_fields           => block_layout_get_num_fields     
     procedure, non_overridable :: get_num_blocks           => block_layout_get_num_blocks
     procedure, non_overridable :: get_block_id             => block_layout_get_block_id
     procedure, non_overridable :: get_block_num_dofs       => block_layout_get_block_num_dofs
     procedure, non_overridable :: get_field_coupling       => block_layout_get_field_coupling
     procedure, non_overridable :: get_field_id_to_block_id => block_layout_get_field_id_to_block_id
     procedure, non_overridable :: get_num_dofs_x_block     => block_layout_get_num_dofs_x_block
     procedure, non_overridable :: fields_coupled           => block_layout_fields_coupled
     procedure, non_overridable :: blocks_coupled           => block_layout_blocks_coupled
     procedure, non_overridable :: set_block_num_dofs       => block_layout_set_block_num_dofs
     procedure, non_overridable :: add_to_block_num_dofs    => block_layout_add_to_block_num_dofs 
     procedure, non_overridable :: equal                    => block_layout_equal
     generic :: operator(==)                                => equal
     
  end type block_layout_t
  
  public :: block_layout_t 
   
contains

  subroutine block_layout_create (this, num_fields, field_id_to_block_id, field_coupling)
    implicit none
    class(block_layout_t)          , intent(inout) :: this
    integer(ip)                    , intent(in)    :: num_fields
    integer(ip)          , optional, intent(in)    :: field_id_to_block_id(num_fields)
    logical              , optional, intent(in)    :: field_coupling(num_fields,num_fields)
    
    integer(ip) :: ifield, jfield
    
    assert ( num_fields >= 1 )
    call this%free()
    this%num_fields = num_fields
    
    if (present(field_id_to_block_id)) then
      assert ( all(field_id_to_block_id >= 1) )
      this%num_blocks = maxval(field_id_to_block_id)
    else
      this%num_blocks = 1
    end if
    
    call memalloc ( this%num_fields, this%field_id_to_block_id, __FILE__, __LINE__ )
    call memalloc ( this%num_blocks, this%num_dofs_x_block, __FILE__, __LINE__ )
    
    if (present(field_id_to_block_id)) then
      this%field_id_to_block_id = field_id_to_block_id
    else
      this%field_id_to_block_id = 1   
    end if
    
    call this%clear_num_dofs_x_block()
    
    call memalloc ( this%num_fields, this%num_fields, this%field_coupling, __FILE__, __LINE__ )
    
    if ( present(field_coupling) ) then
      this%field_coupling = field_coupling
    else
      this%field_coupling = .true.
    end if
    
    call memalloc ( this%num_blocks, this%num_blocks, this%block_coupling, __FILE__, __LINE__ )
    this%block_coupling = .false.
    do jfield=1,this%num_fields
     do ifield=1,this%num_fields
        this%block_coupling(this%field_id_to_block_id(ifield), this%field_id_to_block_id(jfield)) = & 
          this%block_coupling(this%field_id_to_block_id(ifield), this%field_id_to_block_id(jfield)) .or. this%field_coupling(ifield,jfield)
     end do
    end do
  end subroutine block_layout_create
  
  subroutine block_layout_free ( this ) 
    implicit none
    class(block_layout_t), intent(inout) :: this
    this%num_fields = -1
    this%num_blocks = -1 
    if (allocated(this%field_id_to_block_id)) call memfree ( this%field_id_to_block_id, __FILE__, __LINE__ )
    if (allocated(this%num_dofs_x_block)) call memfree ( this%num_dofs_x_block, __FILE__, __LINE__ )
    if (allocated(this%field_coupling)) call memfree ( this%field_coupling, __FILE__, __LINE__ )
    if (allocated(this%block_coupling)) call memfree ( this%block_coupling, __FILE__, __LINE__ )
  end subroutine block_layout_free
  
  subroutine block_layout_clear_num_dofs_x_block ( this ) 
    implicit none
    class(block_layout_t), intent(inout) :: this
    this%num_dofs_x_block = 0
  end subroutine block_layout_clear_num_dofs_x_block

  function block_layout_get_num_fields ( this )
    implicit none
    class(block_layout_t), intent(in) :: this
    integer(ip) :: block_layout_get_num_fields
    block_layout_get_num_fields = this%num_fields
  end function block_layout_get_num_fields
  
  function block_layout_get_num_blocks ( this )
    implicit none
    class(block_layout_t), intent(in) :: this
    integer(ip) :: block_layout_get_num_blocks
    block_layout_get_num_blocks = this%num_blocks
  end function block_layout_get_num_blocks
  
  function block_layout_get_block_id ( this, field_id ) 
    implicit none
    class(block_layout_t), intent(in) :: this
    integer(ip)          , intent(in) :: field_id
    integer(ip) :: block_layout_get_block_id
    assert ( field_id >= 1 .and. field_id <= this%num_fields )
    block_layout_get_block_id = this%field_id_to_block_id(field_id)
  end function block_layout_get_block_id 
  
  pure function block_layout_get_block_num_dofs ( this, block_id ) 
    implicit none
    class(block_layout_t), intent(in) :: this
    integer(ip)          , intent(in) :: block_id
    integer(ip) :: block_layout_get_block_num_dofs
    !assert ( block_id >= 1 .and. block_id <= this%num_blocks )
    block_layout_get_block_num_dofs = this%num_dofs_x_block(block_id)
  end function block_layout_get_block_num_dofs
  
  function block_layout_get_field_coupling ( this ) 
    implicit none
    class(block_layout_t), target, intent(in) :: this
    logical, pointer :: block_layout_get_field_coupling(:,:)
    block_layout_get_field_coupling => this%field_coupling
  end function block_layout_get_field_coupling
  
  function block_layout_get_field_id_to_block_id ( this ) 
    implicit none
    class(block_layout_t), target, intent(in) :: this
    integer(ip), pointer :: block_layout_get_field_id_to_block_id(:)
    block_layout_get_field_id_to_block_id => this%field_id_to_block_id
  end function block_layout_get_field_id_to_block_id
  
  function block_layout_get_num_dofs_x_block ( this ) 
    implicit none
    class(block_layout_t), target, intent(in) :: this
    integer(ip), pointer :: block_layout_get_num_dofs_x_block(:)
    block_layout_get_num_dofs_x_block => this%num_dofs_x_block
  end function block_layout_get_num_dofs_x_block
  
  function block_layout_fields_coupled ( this, ifield, jfield )
    implicit none
    class(block_layout_t), intent(in) :: this
    integer(ip)          , intent(in) :: ifield, jfield
    logical :: block_layout_fields_coupled 
    assert ( ifield >= 1 .and. ifield <= this%num_fields )
    assert ( jfield >= 1 .and. jfield <= this%num_fields )
    block_layout_fields_coupled = this%field_coupling(ifield, jfield)
  end function block_layout_fields_coupled 
  
  function block_layout_blocks_coupled ( this, iblock, jblock )
    implicit none
    class(block_layout_t), intent(in) :: this
    integer(ip)          , intent(in) :: iblock, jblock
    logical :: block_layout_blocks_coupled 
    assert ( iblock >= 1 .and. iblock <= this%num_blocks )
    assert ( jblock >= 1 .and. jblock <= this%num_blocks )
    block_layout_blocks_coupled = this%block_coupling(iblock, jblock)
  end function block_layout_blocks_coupled
  
  subroutine block_layout_set_block_num_dofs ( this, block_id, num_dofs ) 
    implicit none
    class(block_layout_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: block_id
    integer(ip)          , intent(in)    :: num_dofs
    assert ( block_id >= 1 .and. block_id <= this%num_blocks )
    this%num_dofs_x_block(block_id) = num_dofs
  end subroutine block_layout_set_block_num_dofs
  
  subroutine block_layout_add_to_block_num_dofs ( this, block_id, num_dofs ) 
    implicit none
    class(block_layout_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: block_id
    integer(ip)          , intent(in)    :: num_dofs
    assert ( block_id >= 1 .and. block_id <= this%num_blocks )
    this%num_dofs_x_block(block_id) = this%num_dofs_x_block(block_id) + num_dofs
  end subroutine block_layout_add_to_block_num_dofs
  
  function block_layout_equal ( block_layout_op1, block_layout_op2 )
    implicit none
    class(block_layout_t), intent(in) :: block_layout_op1
    class(block_layout_t), intent(in) :: block_layout_op2
    logical :: block_layout_equal
    
    block_layout_equal = (block_layout_op1%num_blocks == block_layout_op2%num_blocks)
    if (.not. block_layout_equal) return
    
    block_layout_equal = (block_layout_op1%num_fields == block_layout_op2%num_fields)  
    if (.not. block_layout_equal) return
    
    block_layout_equal = (allocated(block_layout_op1%field_id_to_block_id) .and. allocated(block_layout_op2%field_id_to_block_id))  
    if (.not. block_layout_equal) return
    
    block_layout_equal = (size(block_layout_op1%field_id_to_block_id) == size(block_layout_op2%field_id_to_block_id))  
    if (.not. block_layout_equal) return
    
    block_layout_equal = (allocated(block_layout_op1%num_dofs_x_block) .and. allocated(block_layout_op2%num_dofs_x_block))  
    if (.not. block_layout_equal) return
    
    block_layout_equal = (size(block_layout_op1%num_dofs_x_block) == size(block_layout_op2%num_dofs_x_block))  
    if (.not. block_layout_equal) return
    
  end function block_layout_equal 
    
end module block_layout_names
