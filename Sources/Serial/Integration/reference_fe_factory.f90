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
module reference_fe_factory_names
  use reference_fe_names
  use types_names
  implicit none
# include "debug.i90"
  private

  public :: make_reference_fe

contains

  function make_reference_fe ( topology, fe_type, number_dimensions, order, field_type, continuity, enable_face_integration )
    implicit none 
    character(*)          , intent(in) :: topology, fe_type
    integer(ip)           , intent(in) :: number_dimensions, order
    character(*)          , intent(in) :: field_type
    logical               , intent(in) :: continuity
    logical , optional    , intent(in) :: enable_face_integration
    
    type(p_reference_fe_t)             :: make_reference_fe
    
    assert ( topology == topology_quad )
    assert ( fe_type  == fe_type_lagrangian )
    
    if ( topology == topology_quad ) then
       if ( fe_type == fe_type_lagrangian ) then
          allocate ( quad_lagrangian_reference_fe_t :: make_reference_fe%p )
       end if
    end if
    call make_reference_fe%p%create( number_dimensions, order, field_type, continuity, enable_face_integration )
  end function make_reference_fe

end module reference_fe_factory_names


