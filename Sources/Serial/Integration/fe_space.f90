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
module fe_space_names
  use types_names
  
  ! Abstract modules
  use assembler_names
  use matrix_array_assembler_names
  use dof_descriptor_names
  use problem_names
  
  implicit none
# include "debug.i90"
  private

  type, abstract :: fe_space_t
     ! A reference to an instance of type(dof_descriptor_t) has
     ! to be available in all subclasses of class(fe_space_t)
     type(dof_descriptor_t), pointer :: dof_descriptor
   contains
     procedure :: set_dof_descriptor => fe_space_set_dof_descriptor
     procedure (volume_integral_interface)                      , deferred :: volume_integral
     procedure (create_matrix_array_assembler_interface)        , deferred :: create_matrix_array_assembler
     procedure (symbolic_setup_matrix_array_assembler_interface), deferred :: symbolic_setup_matrix_array_assembler
  end type fe_space_t
	 
  abstract interface
    ! Selects the dynamic type of class(matrix_array_assembler_t), and that of
    ! the class(matrix_t) and class(array_t) data members that it contains.
    ! It also creates these latter two data members. This subroutine follows
    ! the FACTORY METHOD design pattern.
     function create_matrix_array_assembler_interface(this,& 
          diagonal_blocks_symmetric_storage,&
          diagonal_blocks_symmetric,&
          diagonal_blocks_sign)
       import :: fe_space_t, matrix_array_assembler_t, ip
       implicit none
       class(fe_space_t)              , intent(in) :: this
       logical                        , intent(in) :: diagonal_blocks_symmetric_storage(:)
       logical                        , intent(in) :: diagonal_blocks_symmetric(:)
       integer(ip)                    , intent(in) :: diagonal_blocks_sign(:)
       class(matrix_array_assembler_t), pointer    :: create_matrix_array_assembler_interface
     end function create_matrix_array_assembler_interface
	 
     ! Set-ups symbolically the class(matrix_t) data member that it contains.
     ! Essentially, this amounts for computing the coupling among DoFs in the 
     ! sparsity pattern of class(matrix_t).
     subroutine symbolic_setup_matrix_array_assembler_interface(this,matrix_array_assembler)
       import :: fe_space_t, matrix_array_assembler_t
       implicit none
       class(fe_space_t)              , intent(in)    :: this
       class(matrix_array_assembler_t), intent(inout) :: matrix_array_assembler
     end subroutine symbolic_setup_matrix_array_assembler_interface
	
     subroutine volume_integral_interface(this,approximations,assembler)
       import :: fe_space_t, p_discrete_integration_t,assembler_t
       implicit none
       class(fe_space_t)              , intent(in)    :: this
       class(p_discrete_integration_t), intent(in)    :: approximations(:) 
       class(assembler_t)             , intent(inout) :: assembler
     end subroutine volume_integral_interface
  end interface
  
  ! Data types
  public :: fe_space_t

contains
  subroutine fe_space_set_dof_descriptor(this,dof_descriptor)
    implicit none
    class(fe_space_t), intent(inout) :: this
    type(dof_descriptor_t), target :: dof_descriptor
    this%dof_descriptor => dof_descriptor
  end subroutine fe_space_set_dof_descriptor
end module fe_space_names
