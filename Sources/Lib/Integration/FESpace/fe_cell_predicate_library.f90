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
module fe_cell_predicate_library_names
  use types_names
  use fe_space_names
  
# include "debug.i90"
  implicit none
  private

  type, extends(fe_cell_predicate_t) :: fe_cell_set_id_predicate_t
     private
     integer(ip)    ::    cell_set_id = -1
   contains
     procedure                   ::  visit_fe_cell
     procedure, non_overridable  ::  set_cell_set_id
     procedure, non_overridable  ::  get_cell_set_id
     procedure, non_overridable  ::  free
  end type fe_cell_set_id_predicate_t
  
  public :: fe_cell_set_id_predicate_t
 
contains
  
  ! ========================================================================
  function get_cell_set_id(this)
    implicit none
    class(fe_cell_set_id_predicate_t),  intent(in) :: this
    integer(ip)  :: get_cell_set_id
    get_cell_set_id = this%cell_set_id
  end function get_cell_set_id
  ! ========================================================================
  subroutine set_cell_set_id(this, cell_set_id)
    implicit none
    class(fe_cell_set_id_predicate_t),  intent(inout) :: this
    integer(ip)                      ,  intent(in)    :: cell_set_id
    this%cell_set_id = cell_set_id
  end subroutine set_cell_set_id
  ! ========================================================================
  function visit_fe_cell(this, fe)
    implicit none
    class(fe_cell_set_id_predicate_t),  intent(inout) :: this
    class(fe_cell_iterator_t),          intent(in)    :: fe  
    logical   :: visit_fe_cell
    visit_fe_cell = .false.
    if (fe%get_set_id() == this%cell_set_id) then
       visit_fe_cell = .true.
    end if
  end function visit_fe_cell
  ! ========================================================================
  subroutine free(this)
    implicit none
    class(fe_cell_set_id_predicate_t),  intent(inout) :: this
    this%cell_set_id = -1
  end subroutine free
  
end module fe_cell_predicate_library_names
