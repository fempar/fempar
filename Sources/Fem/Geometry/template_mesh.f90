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
module template_mesh_names
  use types_names
  use migratory_element_names
  use template_element_names
  
  implicit none
  private
  
  type :: template_mesh_t
     integer(ip)                           :: num_elems
     class(migratory_element), allocatable :: mig_elems(:) ! Migratory elements list_t
     type(template_element_t)  ,     pointer :: tmp_elems(:) ! Template elements list_t
   contains
     procedure :: allocate => template_mesh_allocate
     procedure :: free     => template_mesh_free 
     procedure :: print    => template_mesh_print
  end type template_mesh_t
  
  public :: template_mesh_t

contains
  
  subroutine template_mesh_allocate (my, size)
    implicit none
    class(template_mesh_t), intent(out),target :: my
    integer(ip)         , intent(in)  :: size
    
    my%num_elems = size
    allocate( template_element_t :: my%mig_elems(size) )
    select type( this => my%mig_elems )
    type is(template_element_t)
       my%tmp_elems => this
    end select
    
  end subroutine template_mesh_allocate
  
  subroutine template_mesh_free (my)
    implicit none
    class(template_mesh_t), intent(inout) :: my
   
    my%num_elems = -1
    deallocate(my%mig_elems)
    nullify(my%tmp_elems)
    
  end subroutine template_mesh_free

  subroutine template_mesh_print (my)
    implicit none
    class(template_mesh_t), intent(in) :: my
   
    ! Locals
    integer(ip) :: i 

    do i=1, my%num_elems
       write(*,*) 'LID', i, 'GID', my%tmp_elems(i)%globalID, 'area', my%tmp_elems(i)%area
    end do
    
  end subroutine template_mesh_print

end module template_mesh_names
