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
program abstract_reference_element
  implicit none
  type(object) :: object
  !call XXX()
  
  integer(ip), parameter :: DIM = number_space_dimensions
  
  
    
contains
  type object_tree_t
     private
     integer(ip) :: topology
     integer(ip), allocatable :: objects_array(:)     
     integer(ip), allocatable :: ijk_to_index(:)
     integer(ip) :: current = 1    
     contains
     procedure :: create_object_tree
     procedure, private :: fill_tree
     procedure :: create_children_iterator
  end type object_tree_t
  ! Types
  public :: object_tree_t
  
  type children_iterator_t
    private 
    type(object_tree_t), pointer :: object_tree
    integer(ip)                  :: parent
    integer(ip)                  :: i ! Choose a better name for i,j
    integer(ip)                  :: j
  contains
    procedure :: create        => children_iterator_create
    procedure :: init          => children_iterator_init
    procedure :: next          => children_iterator_next
    procedure :: has_finished  => children_iterator_has_finished
    proecudre :: free          => children_iterator_free
  end type children_iterator_t
  
contains




subroutine create_object_tree( this, topology )
   implicit none
   this%topology = topology
   c = 0
   dim = number_of_space_dimensions
   do i = 0,dim
      c = c + 2**(dim-i)*(factorial(dim)/(factorial(i)*factorial(dim-i)))
   end do
   call memalloc ( c, this%objects_array, __FILE__, __LINE__ )       ! c = number objects hypercube dim
   call memalloc ( this%ijk_to_index, ISHFT(ISHFT(1_ip,dim)-1,dim), -1, __FILE__, __LINE__  ) ! (2**dim-1)*2**dim max id object
   
   this%current = 0   
   call this%fill_tree(ISHFT(ISHFT(1_ip,dim)-1,dim)) ! Root object (volume)
   
   call memrealloc ( this%current, this%objects_array, __FILE__, __LINE__ ) 
   call sort( this%current, this%objects_array )
   do i =1,this%current
      ijk_to_index(objects_array(i)) = i
   end do
   this%current = 1
end subroutine

recursive subroutine fill_tree( this, root )
   implicit none
   
   type(children_iterator_t) :: children_iterator
   
   
   if ( this%ijk_to_index(root) /= -1 ) then 
     children_iterator = this%create_children_iterator(root)
     this%current = this%current + 1
     this%objects_array(this%current) = root
     this%ijk_to_index(root) = 1
     do while (.not. children_iterator%has_finished() )
       son = children_iterator%current()
       call this%fill_tree(son)
       call children_iterator%next()
     end do
     !do i = 1, DIM
     !   if ( this%is_component_free( root, i) ) then
     !      do j = 0,1
     !         son = this%generate_son( root, i, j)
     !         if ( this%is_admisible(son) ) then
     !            call this%fill_tree( son )
     !         end if
     !      end do
     !   end if 
     !end do   
     
   end if
   
 end subroutine fill_tree
 
 subroutine is_component_free
    
 subroutine generate_son
 
 subroutine is_admisible
 
 subroutine object_dimension
 
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
  call this%free()
  this%object_tree => object_tree
  this%parent = parent
  call this%init()
end subroutine children_iterator_create

subroutine children_iterator_init ( this)
  implicit none
  class(children_iterator_t), intent(inout) :: this
  this%i = 1
  this%j = 0
  if ( .not. is admisible(i,j) ) then
    call this%children_iterator_next()
  end if  
end subroutine children_iterator_init

recursive subroutine children_iterator_next ( this)
  implicit none
  class(children_iterator_t), intent(inout) :: this
  if ( this%has_finished() ) return 
  if ( this%j == 1 ) then
     this%i = this%i + 1
     this%j = 0
  else
     this%j = 1
  end if
  if ( .not. is admisible(i,j) ) then
      call this%children_iterator_next()
  end if  
end subroutine children_iterator_next  
  
function children_iterator_has_finished ( this)
  implicit none
  class(children_iterator_t), intent(in) :: this
  logical :: children_iterator_has_finished
  children_iterator_has_finished = ( this%i > DIM )
end function children_iterator_has_finished

function children_iterator_current ( this)
  implicit none
  class(children_iterator_t), intent(in) :: this
  integer(ip) :: children_iterator_current
  assert ( .not. this%has_finished() )
  children_iterator_current = F(this%parent,i,j)
end function children_iterator_current


! Creates the son, stores a pointer to the son, and provides the son
              !new_object%free_component(i) = 0 ! Fix component
              !new_object%coordinates(i) = j
              !new_object%object_dimension = this%object_dimension-1
              !new_object%
  






  
  
  
  
contains
  type object_tree
     private
     type(object_t)  :: object_array
     type(has_table) :: bit_to_index_table
     
     integer(ip)     :: object_accesor = 1
     contains
     procedure :: create_object_tree
     procedure :: allocate_object_array
  
  type object_t 
     private 
     ! Data
     ! ...
     integer(ip) :: object_dimension
     integer(ip) :: topology(number_dimensions)
     integer(ip) :: free_components(number_dimensions)
     integer(ip) :: coordinates(number_dimensions)
          
     type(p_object_t), allocatable :: sons()
     
     integer(ip)     :: son_accesor = 1
          
     contains
     ! TBP
     procedure, private :: numbering_all_sub_objects
     procedure, private :: list_sub_objects_of_a_given_dimension
     procedure, private :: apply_displacement_object
     procedure, private :: check_admisible_displacement
  end type object_t
  ! Types
  public :: object_t
contains

  subroutine set_object(this, object_dimension, topology, object_type, coordinates)
     implicit none
     integer(ip) :: object_dimension
     real(rp) :: topology(number_dimensions)
     real(rp) :: object_type(number_dimensions)
     real(rp) :: coordinates(number_dimensions)
     
     this%object_dimension = object_dimension
     this%topology = topology
     this%coordinates = coordinates
     
     call this%allocate_sons()
  end subroutine
    
  subroutine allocate_sons( this )
     implicit none
     allocate ( this%sons, 2**(this%number_dimensions-this%object_dimension), __THIS__, __LINE__ )
  end subroutine
  
  subroutine allocate_object_array( this )
     implicit none
     c = 0
     do i = 0,space_dimension
        c = c + 2**(this%number_dimensions-this%object_dimension)
     end do
     allocate ( this%sons, c, __THIS__, __LINE__ ) ! Maximum number (hexahedron case)
  end subroutine
  
  subroutine start_object_tree(this, object_dimension, topology, object_type)
  
  current_object => this%current_object()
  call create_object( current_object, space_dimension, topology, object_type, coordinates)
  
  end start_object_tree
  
  subroutine create_object_tree( this )
  
  current_object => this%get_current_object()
  do i = 1, space_dimension
     if ( current_object%free_components(i) ) then
        do j = 0,1
           if ( this%is_admissible_son(i,j) ) then
              this%move_object_head()
              son_object => this%get_current_object()
              call current_object%create_son( son_object ) ! Creates the son, stores a pointer to the son, and provides the son
              !new_object%free_component(i) = 0 ! Fix component
              !new_object%coordinates(i) = j
              !new_object%object_dimension = this%object_dimension-1
              !new_object%
              call create_object_tree( this )
              son_object%free()
           end if
        end do
     end if
  end do
  current_object%free()
  
  end
  
  
   

end program abstract_reference_element
