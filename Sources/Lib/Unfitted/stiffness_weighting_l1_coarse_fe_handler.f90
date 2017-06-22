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

module stiffness_weighting_l1_coarse_fe_handler_names
  use fe_space_names
  use par_sparse_matrix_names
  use FPL
  use types_names
  use matrix_names
  use base_sparse_matrix_names
  use par_scalar_array_names
  use environment_names
  use dof_import_names
  use serial_scalar_array_names
  use fe_affine_operator_names
  use base_static_triangulation_names
  use cell_import_names

  implicit none
# include "debug.i90"
  private

!========================================================================================
  type, extends(standard_l1_coarse_fe_handler_t) :: stiffness_weighting_l1_coarse_fe_handler_t
    private

    class(par_sparse_matrix_t), pointer :: matrix => null()

  contains

    ! Added TBPs
    generic   :: create                   => stiffness_l1_create
    procedure, private :: stiffness_l1_create
    procedure :: free                     => stiffness_l1_free

    ! Overwritten TBPs
    procedure :: setup_weighting_operator => stiffness_l1_setup_weighting_operator

  end type stiffness_weighting_l1_coarse_fe_handler_t

  public :: stiffness_weighting_l1_coarse_fe_handler_t

contains

!========================================================================================
subroutine stiffness_l1_create(this, fe_affine_operator)

  implicit none
  class(stiffness_weighting_l1_coarse_fe_handler_t), intent(inout) :: this
  class(fe_affine_operator_t), target,    intent(in)    :: fe_affine_operator

  class(matrix_t),            pointer :: matrix

  call this%free()

  matrix => fe_affine_operator%get_matrix()
  select type (matrix)
    class is (par_sparse_matrix_t)
      this%matrix => matrix
    class default
      check(.false.)
  end select

end subroutine stiffness_l1_create

!========================================================================================
subroutine stiffness_l1_free(this)
  implicit none
  class(stiffness_weighting_l1_coarse_fe_handler_t), intent(inout) :: this
  this%matrix         => null()
end subroutine stiffness_l1_free

!========================================================================================
subroutine stiffness_l1_setup_weighting_operator(this,field_id,par_fe_space,parameter_list,weighting_operator)
  implicit none
  class(stiffness_weighting_l1_coarse_fe_handler_t), intent(in)    :: this
    integer(ip)                         , intent(in)    :: field_id
  type(par_fe_space_t)                  , intent(in)    :: par_fe_space
  type(parameterlist_t)                 , intent(in)    :: parameter_list
  real(rp), allocatable                 , intent(inout) :: weighting_operator(:)

  integer(ip)                          :: block_id
  integer(ip), pointer                 :: field_to_block(:)
  integer(ip)                          :: num_dofs
  type(par_scalar_array_t)             :: par_array
  type(environment_t), pointer         :: p_env
  type(dof_import_t),  pointer         :: dof_import
  type(serial_scalar_array_t), pointer :: serial_array
  real(rp), pointer                    :: assembled_diag(:)
  real(rp), allocatable                :: sub_assembled_diag(:)
  integer(ip)                          :: istat


  ! We assume a single field for the moment
  assert(field_id == 1)
  assert(par_fe_space%get_number_fields() == 1)

  ! Clean up
  if (allocated(weighting_operator) ) then
    call memfree ( weighting_operator, __FILE__, __LINE__ )
  end if

  ! Allocate the weighting
  field_to_block => par_fe_space%get_field_blocks()
  block_id = field_to_block(field_id)
  num_dofs = par_fe_space%get_block_number_dofs(block_id)
  call memalloc(num_dofs,weighting_operator,__FILE__,__LINE__)

  ! Get the sub-assembled diagonal
  massert(associated(this%matrix),'The matrix is not associated. Have you called `stiffness_weighting_l1_coarse_fe_handler_t%create` ? ')
  call this%matrix%extract_diagonal(sub_assembled_diag)

  ! Communicate to compute the fully assembled diagonal
  p_env => par_fe_space%get_environment()
  dof_import => par_fe_space%get_block_dof_import(block_id)
  call par_array%create_and_allocate(p_env, dof_import)
  serial_array   => par_array%get_serial_scalar_array()
  assembled_diag => serial_array%get_entries()
  assembled_diag(:) = sub_assembled_diag(:)
  call par_array%comm()

  ! Compute the weighting
  weighting_operator(:) = sub_assembled_diag(:)/assembled_diag(:)

  ! Clean up
  deallocate(sub_assembled_diag,stat=istat); check(istat == 0)
  call par_array%free()

end subroutine stiffness_l1_setup_weighting_operator

end module stiffness_weighting_l1_coarse_fe_handler_names
!***************************************************************************************************
