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
module fe_affine_operator_names
  use types_names
  use memor_names
  use problem_names
  
  ! Abstract types
  use dof_descriptor_names
  use fe_space_names
  use operator_names
  use vector_names
  use matrix_array_assembler_names
  use array_names
  use matrix_names

  implicit none
# include "debug.i90"

  private
  

  type, extends(operator_t):: fe_affine_operator_t
    private
    class(fe_space_t)              , pointer :: fe_space               => NULL()
	class(p_discrete_integration_t), pointer :: approximations(:)      => NULL()
	class(matrix_array_assembler_t), pointer :: matrix_array_assembler => NULL()
  contains
    procedure :: create          => fe_affine_operator_create
	procedure :: symbolic_setup  => fe_affine_operator_symbolic_setup
	procedure :: numerical_setup => fe_affine_operator_numerical_setup
	procedure :: apply           => fe_affine_operator_apply
	procedure :: apply_fun       => fe_affine_operator_apply_fun
	procedure :: get_matrix      => fe_affine_operator_get_matrix
	procedure :: get_array       => fe_affine_operator_get_array
	procedure :: free_in_stages  => fe_affine_operator_free_in_stages
	procedure :: free            => fe_affine_operator_free
  end type fe_affine_operator_t

  ! Types
  public :: fe_affine_operator_t

contains
  subroutine fe_affine_operator_create (this,&
										diagonal_blocks_symmetric_storage,&
										diagonal_blocks_symmetric,&
										diagonal_blocks_sign,&
										fe_space,&
										approximations)
    implicit none
	class(fe_affine_operator_t)            , intent(out) :: this
	class(fe_space_t)              , target, intent(in)  :: fe_space
    logical                                , intent(in)  :: diagonal_blocks_symmetric_storage(fe_space%dof_descriptor%nblocks)
    logical                                , intent(in)  :: diagonal_blocks_symmetric(fe_space%dof_descriptor%nblocks)
    integer(ip)                            , intent(in)  :: diagonal_blocks_sign(fe_space%dof_descriptor%nblocks)
	class(p_discrete_integration_t), target, intent(in)  :: approximations(:)
	
	this%fe_space               => fe_space
    this%approximations         => approximations
    this%matrix_array_assembler => fe_space%create_matrix_array_assembler(diagonal_blocks_symmetric_storage, &
																		  diagonal_blocks_symmetric, &
																		  diagonal_blocks_sign)					
  end subroutine fe_affine_operator_create
  
  subroutine fe_affine_operator_symbolic_setup (this)
    implicit none
	class(fe_affine_operator_t), intent(inout) :: this
	assert ( associated(this%fe_space) )
	assert ( associated(this%matrix_array_assembler) )
	call this%fe_space%symbolic_setup_matrix_array_assembler(this%matrix_array_assembler)
  end subroutine fe_affine_operator_symbolic_setup
  
  subroutine fe_affine_operator_numerical_setup (this)
    implicit none
	class(fe_affine_operator_t), intent(inout) :: this
	assert ( associated(this%fe_space) )
	assert ( associated(this%matrix_array_assembler) )
	call this%matrix_array_assembler%allocate()
	call this%fe_space%volume_integral(this%approximations,this%matrix_array_assembler)
  end subroutine fe_affine_operator_numerical_setup
 
  subroutine fe_affine_operator_free_in_stages(this,action)
    implicit none
	class(fe_affine_operator_t), intent(inout) :: this
	integer(ip)                , intent(in)    :: action
	integer(ip)                                :: istat
	assert ( associated(this%fe_space) )
	assert ( associated(this%matrix_array_assembler) )
	call this%matrix_array_assembler%free_in_stages(action)
	if ( action == free_clean ) then
	   deallocate ( this%matrix_array_assembler, stat=istat )
	   check(istat==0)
	   nullify(this%matrix_array_assembler)
	   nullify(this%fe_space)
       nullify(this%approximations)
	end if
  end subroutine fe_affine_operator_free_in_stages
  
  subroutine fe_affine_operator_free(this)
    implicit none
	class(fe_affine_operator_t), intent(inout) :: this
    call this%free_in_stages(free_values)
	call this%free_in_stages(free_struct)
	call this%free_in_stages(free_clean)
  end subroutine fe_affine_operator_free
  
  function fe_affine_operator_get_matrix(this)
    implicit none
	class(fe_affine_operator_t), target, intent(in) :: this
    class(matrix_t), pointer :: fe_affine_operator_get_matrix
	fe_affine_operator_get_matrix => this%matrix_array_assembler%get_matrix()
  end function fe_affine_operator_get_matrix
  
  function fe_affine_operator_get_array(this)
    implicit none
	class(fe_affine_operator_t), target, intent(in) :: this
	class(array_t), pointer :: fe_affine_operator_get_array
    fe_affine_operator_get_array => this%matrix_array_assembler%get_array()
  end function fe_affine_operator_get_array
  
  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine fe_affine_operator_apply(op,x,y) 
    implicit none
    class(fe_affine_operator_t), intent(in)    :: op
    class(vector_t) , intent(in)    :: x
    class(vector_t) , intent(inout) :: y 
	class(matrix_t) , pointer       :: matrix
	class(array_t)  , pointer       :: array

    call x%GuardTemp()
    matrix => op%matrix_array_assembler%get_matrix()
	call matrix%apply(x,y)
    call x%CleanTemp()
  end subroutine fe_affine_operator_apply
  
  ! op%apply(x)
  ! Allocates room for (temporary) y
  function fe_affine_operator_apply_fun(op,x) result(y)
    implicit none
    class(fe_affine_operator_t), intent(in)  :: op
    class(vector_t) , intent(in)  :: x
    class(vector_t) , allocatable :: y 
	class(matrix_t), pointer :: matrix

	call x%GuardTemp()
	matrix => op%matrix_array_assembler%get_matrix()
	allocate(y, mold=x); call y%default_initialization()
    y = matrix%apply_fun(x) 
	call x%CleanTemp()
	call y%SetTemp()
  end function fe_affine_operator_apply_fun
  
end module fe_affine_operator_names
