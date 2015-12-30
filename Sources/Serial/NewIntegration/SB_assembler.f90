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
module SB_assembler_names
  use types_names
  use dof_descriptor_names

  implicit none
# include "debug.i90"
  private

  type, abstract :: SB_assembler_t
    contains
	     procedure (assembly_interface)        , deferred :: assembly
      procedure (compress_storage_interface), deferred :: compress_storage
  end type
  
  abstract interface
     subroutine assembly_interface( this, & 
                                    number_fe_spaces, &
                                    number_nodes, &
                                    elem2dof, &
                                    field_blocks, &
                                    field_coupling, &
                                    elmat, &
                                    elvec )
       import :: SB_assembler_t, rp, ip, i1p_t
       implicit none
       class(SB_assembler_t) , intent(inout) :: this
       integer(ip)           , intent(in)    :: number_fe_spaces
       integer(ip)           , intent(in)    :: number_nodes(number_fe_spaces)
       type(i1p_t)           , intent(in)    :: elem2dof(number_fe_spaces)
       integer(ip)           , intent(in)    :: field_blocks(number_fe_spaces)
       logical               , intent(in)    :: field_coupling(number_fe_spaces,number_fe_spaces)
       ! elmat MUST have as many rows/columns as \sum_{i=1}^{number_fe_spaces} number_nodes(i)
       real(rp)              , intent(in)    :: elmat(:,:) 
       ! elvec MUST have as many entries as \sum_{i=1}^{number_fe_spaces} number_nodes(i)
       real(rp)              , intent(in)    :: elvec(:)   
     end subroutine assembly_interface
     
     subroutine compress_storage_interface( this, & 
                                            sparse_matrix_storage_format )
       import :: SB_assembler_t
       implicit none
       class(SB_assembler_t) , intent(inout) :: this
       character(*)          , intent(in)    :: sparse_matrix_storage_format
     end subroutine compress_storage_interface
     
  end interface
	 
  ! Data types
  public :: SB_assembler_t

end module SB_assembler_names
