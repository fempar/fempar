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
module par_assembly_names
  use types_names
  use assembly_names
  use finite_element_names
  use integrable_names
  use dof_descriptor_names
  use par_scalar_matrix_names
  use par_scalar_array_names
  use par_block_matrix_names
  use par_block_array_names
  use par_scalar_names
  
  implicit none
# include "debug.i90"
  private
  
  public :: par_assembly

contains

  subroutine par_assembly(finite_element, dof_descriptor, a ) 
    implicit none
    ! Parameters
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    class(integrable_t)   , intent(inout) :: a

    select type(a)
    class is(par_scalar_matrix_t)
       call assembly_element_par_matrix_mono(finite_element, dof_descriptor, a) 
    class is(par_scalar_array_t)
       call assembly_element_par_vector_mono(finite_element, dof_descriptor, a)
    class is(par_block_matrix_t)
       call assembly_element_par_matrix_block(finite_element, dof_descriptor, a)
    class is(par_block_array_t)
       call assembly_element_par_vector_block(finite_element, dof_descriptor, a)
    class is(par_scalar_t)
       call assembly_element_par_scalar(finite_element, dof_descriptor, a)
    class default
       ! class not yet implemented
       check(.false.)
    end select
  end subroutine par_assembly

  subroutine assembly_element_par_matrix_block(  finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t)  , intent(in)    :: dof_descriptor
    type(finite_element_t)  , intent(in)    :: finite_element
    type(par_block_matrix_t), intent(inout) :: a
    
    integer(ip) :: ivar, iblock, jblock

    type(par_scalar_matrix_t), pointer :: p_matrix

    do iblock = 1, dof_descriptor%nblocks
       do jblock = 1, dof_descriptor%nblocks
          p_matrix => a%get_block(iblock,jblock)
          if ( associated(p_matrix) ) then
             if(p_matrix%p_env%am_i_fine_task()) then
                call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, p_matrix%f_matrix, iblock, jblock )
             end if
          end if 
       end do
    end do

  end subroutine assembly_element_par_matrix_block

  subroutine assembly_element_par_matrix_mono(  finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(par_scalar_matrix_t)    , intent(inout) :: a

    if(a%p_env%am_i_fine_task()) then
       call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, a%f_matrix )
    end if

  end subroutine assembly_element_par_matrix_mono

  subroutine assembly_element_par_vector_block(  finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t)  , intent(in)    :: dof_descriptor
    type(finite_element_t)  , intent(in)    :: finite_element
    type(par_block_array_t), intent(inout) :: a

    integer(ip) :: iblock

    do iblock = 1, dof_descriptor%nblocks
       if(a%blocks(iblock)%p_env%am_i_fine_task()) then
          call element_serial_scalar_array_assembly( dof_descriptor, finite_element, a%blocks(iblock)%f_vector, &
               & iblock )
       end if
    end do

  end subroutine assembly_element_par_vector_block

  subroutine assembly_element_par_vector_mono(  finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(par_scalar_array_t)    , intent(inout) :: a

    if(a%p_env%am_i_fine_task()) then
       call element_serial_scalar_array_assembly( dof_descriptor, finite_element, a%f_vector )
    end if

  end subroutine assembly_element_par_vector_mono

  subroutine assembly_element_par_scalar ( finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(finite_element_t), intent(in)    :: finite_element
    type(par_scalar_t)    , intent(inout) :: a
    
    call a%sum(finite_element%scalar)

  end subroutine assembly_element_par_scalar

end module par_assembly_names
