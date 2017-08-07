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
module hp_adaptive_fe_space_names
  use types_names
  use memor_names
  use base_static_triangulation_names
  use p4est_serial_triangulation_names
  use reference_fe_names
  use fe_space_names
  use fe_function_names
  use environment_names
  use conditions_names
  use std_vector_integer_ip_names
  use std_vector_real_rp_names
  use list_types_names
  
  use matrix_array_assembler_names
  use serial_scalar_array_names
  use vector_names
  use array_names
  use matrix_names
  use sparse_matrix_names
  use block_layout_names
  
  
  implicit none
# include "debug.i90"
  private

  ! TODO (potentially): to move member variables of serial_hp_adaptive_fe_space_t to serial_fe_space_t
  type, extends(serial_fe_space_t) :: serial_hp_adaptive_fe_space_t
     private
     !list_t :: constraints 
     ! Hanging and strong Dirichlet data
     ! ptr_constraint_dofs         : pointer to constraints
     ! l1                          : constraint DOFs dependencies (0 for independent term)
     ! constraint_dofs_coefficients: constraint DoFs coefficients (also independent term)
     ! u_fixed = sum u_dep w_dep + c
     integer(ip)                                 :: number_fixed_dofs = -1
     type(std_vector_integer_ip_t)               :: ptr_constraint_dofs
     type(std_vector_integer_ip_t)               :: constraint_dofs_dependencies
     type(std_vector_real_rp_t)                  :: constraint_dofs_coefficients
     
     type(p4est_serial_triangulation_t), pointer :: p4est_triangulation =>  NULL()
   contains
     procedure          :: create_fe_vef_iterator                                 => serial_hp_adaptive_fe_space_create_fe_vef_iterator
     procedure          :: create_fe_iterator                                     => serial_hp_adaptive_fe_space_create_fe_iterator
     procedure          :: create_fe_face_iterator                                => serial_hp_adaptive_fe_space_create_fe_face_iterator
     
     procedure          :: serial_fe_space_create_same_reference_fes_on_all_cells => shpafs_create_same_reference_fes_on_all_cells 
     procedure          :: serial_fe_space_create_different_between_cells         => shpafs_create_different_between_cells
     procedure          :: free                                                   => serial_hp_adaptive_fe_space_free
     
     procedure          :: fill_dof_info                                          => serial_hp_adaptive_fe_space_fill_dof_info
     procedure, private :: fill_elem2dof_and_count_dofs                           => serial_hp_adaptive_fe_space_fill_elem2dof_and_count_dofs
     
     procedure          :: setup_hanging_node_constraints                         => shpafs_setup_hanging_node_constraints
     procedure          :: transfer_dirichlet_to_constraint_dof_coefficients      => shpafs_transfer_dirichlet_to_constraint_dof_coefficients
     procedure          :: transfer_dirichlet_to_fe_space                         => shpafs_transfer_dirichlet_to_fe_space
     procedure          :: free_ptr_constraint_dofs                               => shpafs_free_ptr_constraint_dofs
     procedure          :: free_constraint_dofs_dependencies                      => shpafs_free_constraint_dofs_dependencies
     procedure          :: free_constraint_dofs_coefficients                      => shpafs_free_constraint_dofs_coefficients
     
     procedure          :: get_number_fixed_dofs                                  => shpafs_get_number_fixed_dofs
     procedure          :: set_up_strong_dirichlet_bcs                            => shpafs_set_up_strong_dirichlet_bcs
     procedure          :: update_fixed_dof_values                                => shpafs_update_fixed_dof_values
     procedure          :: interpolate_dirichlet_values                           => shpafs_interpolate_dirichlet_values
     
     procedure, private :: fill_vef_lids_of_fe_faces                              => shpafs_fill_vef_lids_of_fe_faces
     
     procedure          :: project_ref_fe_id_per_fe                               => shpafs_project_ref_fe_id_per_fe
     procedure          :: project_fe_integration_arrays                          => shpafs_project_fe_integration_arrays
     procedure          :: project_fe_face_integration_arrays                     => shpafs_project_fe_face_integration_arrays
     procedure          :: refine_and_coarsen                                     => serial_hp_adaptive_fe_space_refine_and_coarsen
 end type serial_hp_adaptive_fe_space_t  
 
 type, extends(fe_iterator_t) :: hp_adaptive_fe_iterator_t
   private 
   type(serial_hp_adaptive_fe_space_t), pointer :: hp_adaptive_fe_space => NULL()
 contains
   procedure          :: create                     => hp_adaptive_fe_iterator_create
   procedure          :: free                       => hp_adaptive_fe_iterator_free
   procedure          :: assemble                   => hp_adaptive_fe_iterator_assemble
   procedure, private :: recursive_matrix_assembly  => hp_adaptive_fe_iterator_recursive_matrix_assembly
   procedure, private :: recursive_vector_assembly  => hp_adaptive_fe_iterator_recursive_vector_assembly
 end type hp_adaptive_fe_iterator_t
 
 type, extends(fe_face_iterator_t) :: hp_adaptive_fe_face_iterator_t
   private 
   type(serial_hp_adaptive_fe_space_t), pointer :: hp_adaptive_fe_space => NULL()
 contains
   procedure          :: create                     => hp_adaptive_fe_face_iterator_create
   procedure          :: free                       => hp_adaptive_fe_face_iterator_free
   procedure          :: assemble                   => hp_adaptive_fe_face_iterator_assemble
   procedure, private :: recursive_matrix_assembly  => hp_adaptive_fe_face_iterator_recursive_matrix_assembly
   procedure, private :: recursive_vector_assembly  => hp_adaptive_fe_face_iterator_recursive_vector_assembly
 end type hp_adaptive_fe_face_iterator_t
 
 public :: serial_hp_adaptive_fe_space_t
 
contains


subroutine hp_adaptive_fe_iterator_create ( this, fe_space )
  implicit none
  class(hp_adaptive_fe_iterator_t), intent(inout) :: this
  class(serial_fe_space_t), target, intent(in)    :: fe_space
  call this%free()
  call this%set_fe_space(fe_space)
  select type(fe_space)
  class is (serial_hp_adaptive_fe_space_t)
    this%hp_adaptive_fe_space => fe_space
  class default
    assert(.false.)
  end select
  call this%create_cell(this%hp_adaptive_fe_space%p4est_triangulation)
end subroutine hp_adaptive_fe_iterator_create

subroutine hp_adaptive_fe_iterator_free (this)
  implicit none
  class(hp_adaptive_fe_iterator_t), intent(inout) :: this
  if ( associated(this%hp_adaptive_fe_space ) ) then
    if ( associated(this%hp_adaptive_fe_space%p4est_triangulation) ) then
      call this%free_cell(this%hp_adaptive_fe_space%p4est_triangulation)
    end if
  end if
  call this%nullify_fe_space()
  nullify(this%hp_adaptive_fe_space)
end subroutine hp_adaptive_fe_iterator_free

!! Assembly of local matrices for hp-adaptivity
subroutine hp_adaptive_fe_iterator_assemble(this,elmat,elvec,matrix_array_assembler)
  implicit none
  class(hp_adaptive_fe_iterator_t), intent(in)    :: this
  real(rp)                        , intent(in)    :: elmat(:,:)
  real(rp)                        , intent(in)    :: elvec(:)
  class(matrix_array_assembler_t) , intent(inout) :: matrix_array_assembler
  
  integer(ip), pointer :: elem2dof(:)
  integer(ip) :: i, j
  
  class(matrix_t), pointer :: matrix
  class(array_t), pointer :: array
  type(sparse_matrix_t), pointer :: sparse_matrix
  type(serial_scalar_array_t), pointer :: serial_scalar_array
  
  
  call this%get_field_elem2dof(1,elem2dof)
  matrix => matrix_array_assembler%get_matrix()
  array  => matrix_array_assembler%get_array()  
  
  select type(matrix)
  class is(sparse_matrix_t)
     sparse_matrix => matrix
  class default
     assert(.false.)
  end select
  
  select type(array)
  class is(serial_scalar_array_t)
     serial_scalar_array => array
  class default
     assert(.false.)
  end select
  
  do j=1, this%get_number_dofs()
     do i=1, this%get_number_dofs() 
        call this%recursive_matrix_assembly( elem2dof(i), elem2dof(j), elmat(i,j), sparse_matrix, serial_scalar_array)
     end do
     call this%recursive_vector_assembly (elem2dof(j), elvec(j), serial_scalar_array)
  end do 
end subroutine hp_adaptive_fe_iterator_assemble

recursive subroutine hp_adaptive_fe_iterator_recursive_matrix_assembly(this, i, j, a_ij, matrix, array)
  implicit none
  class(hp_adaptive_fe_iterator_t), intent(in)    :: this
  integer(ip)                     , intent(in)    :: i
  integer(ip)                     , intent(in)    :: j
  real(rp)                        , intent(in)    :: a_ij
  type(sparse_matrix_t)           , intent(inout) :: matrix
  type(serial_scalar_array_t)     , intent(inout) :: array
  
  integer(ip) :: k, pos, spos, epos
  real(rp) :: weight
 
  if ( .not. this%is_fixed_dof(i)) then    
    if ( .not. this%is_fixed_dof(j) ) then 
      ! Insert a_ij, v_j on matrix_array_assembler position (i,j)  
      call matrix%insert(i,j,a_ij)
    else ! j is a fixed DoF
      if ( this%is_strong_dirichlet_dof(j) ) then
        assert ( this%hp_adaptive_fe_space%constraint_dofs_dependencies%get(this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(j))) == 0 )
        weight  = this%hp_adaptive_fe_space%constraint_dofs_coefficients%get(this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(j)))
        call array%add(i,-a_ij*weight)
      else
        ! Traverse DoFs on which j depends on
        spos = this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(j))
        epos = this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(j)+1)-1
        do pos=spos, epos 
          k       = this%hp_adaptive_fe_space%constraint_dofs_dependencies%get(pos)
          weight  = this%hp_adaptive_fe_space%constraint_dofs_coefficients%get(pos)
          call this%recursive_matrix_assembly( i, k, weight*a_ij, matrix, array)
        end do 
       end if 
      end if
   else ! i is a fixed DoF
     if ( .not. this%is_strong_dirichlet_dof(i) ) then
        ! Traverse DoFs on which i depends on
        spos = this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(i))
        epos = this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(i)+1)-1
        do pos=spos, epos 
          k       = this%hp_adaptive_fe_space%constraint_dofs_dependencies%get(pos)
          weight  = this%hp_adaptive_fe_space%constraint_dofs_coefficients%get(pos)
          call this%recursive_matrix_assembly( k, j, weight*a_ij, matrix, array)
        end do
      end if
  end if
end subroutine hp_adaptive_fe_iterator_recursive_matrix_assembly

recursive subroutine hp_adaptive_fe_iterator_recursive_vector_assembly(this, i, v_i, array)
  implicit none
  class(hp_adaptive_fe_iterator_t), intent(in)    :: this
  integer(ip)                     , intent(in)    :: i
  real(rp)                        , intent(in)    :: v_i
  type(serial_scalar_array_t)     , intent(inout) :: array
  
  integer(ip) :: k, pos, spos, epos
  real(rp) :: weight
 
  if ( .not. this%is_fixed_dof(i) ) then    
      call array%add(i,v_i)
  else ! i is as fixed DoF
    if ( .not. this%is_strong_dirichlet_dof(i) ) then
        ! Traverse DoFs on which i depends on
        spos = this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(i))
        epos = this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(i)+1)-1
        do pos=spos, epos 
          k       = this%hp_adaptive_fe_space%constraint_dofs_dependencies%get(pos)
          weight  = this%hp_adaptive_fe_space%constraint_dofs_coefficients%get(pos)
          call this%recursive_vector_assembly( k, weight*v_i, array )
        end do
     end if
  end if
end subroutine hp_adaptive_fe_iterator_recursive_vector_assembly


!  
!subroutine fe_iterator_recursive_assembly(this, i, j, a_ij, v_j, matrix, array)
!  implicit none
!  class(fe_iterator_t)           , intent(in)    :: this
!  integer(ip)                    , intent(in)    :: i
!  integer(ip)                    , intent(in)    :: j
!  real(rp)                       , intent(in)    :: a_ij
!  real(rp)                       , intent(in)    :: v_j
!  class(matrix_t)                , intent(inout) :: matrix
!  class(array_t)                 , intent(inout) :: array
!  
!  if (i < 0) then
!      ! Traverse DoFs on which i depends on
!      do pos=this%ptr_constraint(abs(i)), this%ptr_constraint(abs(i)+1)-1
!        k       = this%constraint_dofs_dependencies%get(pos)
!        if (k>0) then ! If current DoF on which i depends is subject to Dirichlet BC's
!          weight  = this%constraint_coefficient%get(pos)
!          call this%recursive_assembly( k, j, weight*a_ij, v_j, matrix_array_assembler )
!        end if  
!      end do
!  else
!    if ( j < 0 ) then
!      ! Traverse DoFs on which j depends on
!      do pos=this%ptr_constraint(abs(j)), this%ptr_constraint(abs(j)+1)-1
!        k       = this%constraint_dofs_dependencies%get(pos)
!        weight  = this%constraint_coefficient%get(pos)
!        if (k>0) then ! If current DoF on which i depends is subject to Dirichlet BC's
!          call this%recursive_assembly( i, k, weight*a_ij, weight*v_j, matrix_array_assembler )
!        else
!          !  Insert v_j-a_ij*weight in position j of RHS array
!          select type(array)
!          type is (serial_scalar_array_t)
!            call array%add_entry(j,v_j-a_ij*weight)
!          type is (par_scalar_array_t)
!            call array%add_entry(j,v_j-a_ij*weight)
!          default
!            assert(.false.)
!          end select
!        end if
!      end do   
!    else
!        ! Insert a_ij, v_j on matrix_array_assembler position (i,j)  
!        select type(matrix)
!          type is (sparse_matrix_t)
!            call matrix%insert(i,j,a_ij)
!          type is (par_sparse_matrix_t)
!            call matrix%insert(i,j,a_ij)
!          default 
!            assert(.false.)
!        end select
!        
!        select type(array)
!          type is (serial_scalar_array_t)
!            call array%add_entry(j,v_j)
!          type is (par_scalar_array_t)
!            call array%add_entry(j,v_j)
!          default: 
!            assert(.false.)
!        end select
!    end if
!  end if
!end fe_iterator_recursive_assembly

subroutine hp_adaptive_fe_face_iterator_create ( this, fe_space, vef )
  implicit none
  class(hp_adaptive_fe_face_iterator_t), target, intent(inout) :: this
  class(serial_fe_space_t)             , target, intent(in)    :: fe_space
  class(vef_iterator_t)                        , intent(in)    :: vef
  select type(fe_space)
  class is (serial_hp_adaptive_fe_space_t)
    this%hp_adaptive_fe_space => fe_space
  class default
    assert(.false.)
  end select
  call this%fe_face_iterator_t%create(fe_space,vef)
end subroutine hp_adaptive_fe_face_iterator_create

subroutine hp_adaptive_fe_face_iterator_free ( this )
  implicit none
  class(hp_adaptive_fe_face_iterator_t), intent(inout) :: this
  integer(ip) :: istat
  call this%fe_face_iterator_t%free()
  nullify(this%hp_adaptive_fe_space)
end subroutine hp_adaptive_fe_face_iterator_free

subroutine hp_adaptive_fe_face_iterator_assemble(this,facemat,facevec,matrix_array_assembler)
  implicit none
  class(hp_adaptive_fe_face_iterator_t), intent(in)    :: this
  real(rp)                             , intent(in)    :: facemat(:,:,:,:)
  real(rp)                             , intent(in)    :: facevec(:,:)
  class(matrix_array_assembler_t)      , intent(inout) :: matrix_array_assembler
  
  class(serial_fe_space_t)       , pointer     :: fe_space
  type(hp_adaptive_fe_iterator_t)              :: test_fe, trial_fe
  integer(ip)                                  :: ineigh, jneigh, i, j
  integer(ip)                                  :: test_number_dofs, trial_number_dofs
  integer(ip)                    , pointer     :: test_elem2dof(:), trial_elem2dof(:)
  
  class(matrix_t)            , pointer :: matrix
  class(array_t)             , pointer :: array
  type(sparse_matrix_t)      , pointer :: sparse_matrix
  type(serial_scalar_array_t), pointer :: serial_scalar_array
  
  matrix => matrix_array_assembler%get_matrix()
  array  => matrix_array_assembler%get_array()  
  
  select type(matrix)
  class is(sparse_matrix_t)
    sparse_matrix => matrix
  class default
    assert(.false.)
  end select
  
  select type(array)
  class is(serial_scalar_array_t)
    serial_scalar_array => array
  class default
    assert(.false.)
  end select
  
  fe_space => this%get_fe_space()
  call test_fe%create(fe_space)
  call trial_fe%create(fe_space)
  
  do ineigh = 1,this%get_num_cells_around()
    call this%get_cell_around(ineigh,test_fe)
    test_number_dofs = test_fe%get_number_dofs()
    call test_fe%get_field_elem2dof(1,test_elem2dof)
    do jneigh = 1,this%get_num_cells_around()
      call this%get_cell_around(jneigh,trial_fe)
      trial_number_dofs = trial_fe%get_number_dofs()
      call trial_fe%get_field_elem2dof(1,trial_elem2dof)
      do i = 1,test_number_dofs
         do j = 1,trial_number_dofs
           call this%recursive_matrix_assembly( test_fe,                    &
                                                trial_fe,                   &
                                                test_elem2dof(i),           &
                                                trial_elem2dof(j),          &
                                                facemat(i,j,ineigh,jneigh), &
                                                sparse_matrix,              &
                                                serial_scalar_array )
         end do
      end do
    end do
    do i = 1,test_number_dofs
      call this%recursive_vector_assembly( test_fe,           &
                                           test_elem2dof(i),  &
                                           facevec(i,ineigh), &
                                           serial_scalar_array )
    end do
  end do
  
  call trial_fe%free()
  call test_fe%free()
  
end subroutine hp_adaptive_fe_face_iterator_assemble

recursive subroutine hp_adaptive_fe_face_iterator_recursive_matrix_assembly(this, test_fe, trial_fe, i, j, a_ij, matrix, array)
  implicit none
  class(hp_adaptive_fe_face_iterator_t), intent(in)    :: this
  type(hp_adaptive_fe_iterator_t)      , intent(in)    :: test_fe
  type(hp_adaptive_fe_iterator_t)      , intent(in)    :: trial_fe
  integer(ip)                          , intent(in)    :: i
  integer(ip)                          , intent(in)    :: j
  real(rp)                             , intent(in)    :: a_ij
  type(sparse_matrix_t)                , intent(inout) :: matrix
  type(serial_scalar_array_t)          , intent(inout) :: array
  
  integer(ip) :: k, pos, spos, epos
  real(rp) :: weight
  
  if ( .not. test_fe%is_fixed_dof(i)) then    
    if ( .not. trial_fe%is_fixed_dof(j) ) then 
      ! Insert a_ij, v_j on matrix_array_assembler position (i,j)  
      call matrix%insert(i,j,a_ij)
    else ! j is a fixed DoF
      if ( trial_fe%is_strong_dirichlet_dof(j) ) then
        assert ( trial_fe%hp_adaptive_fe_space%constraint_dofs_dependencies%get(this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(j))) == 0 )
        weight  = trial_fe%hp_adaptive_fe_space%constraint_dofs_coefficients%get(this%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(j)))
        call array%add(i,-a_ij*weight)
      else
        ! Traverse DoFs on which j depends on
        spos = trial_fe%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(j))
        epos = trial_fe%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(j)+1)-1
        do pos=spos, epos 
          k       = trial_fe%hp_adaptive_fe_space%constraint_dofs_dependencies%get(pos)
          weight  = trial_fe%hp_adaptive_fe_space%constraint_dofs_coefficients%get(pos)
          call this%recursive_matrix_assembly( test_fe, trial_fe, i, k, weight*a_ij, matrix, array)
        end do 
       end if 
      end if
   else ! i is a fixed DoF
     if ( .not. test_fe%is_strong_dirichlet_dof(i) ) then
        ! Traverse DoFs on which i depends on
        spos = test_fe%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(i))
        epos = test_fe%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(i)+1)-1
        do pos=spos, epos 
          k       = test_fe%hp_adaptive_fe_space%constraint_dofs_dependencies%get(pos)
          weight  = test_fe%hp_adaptive_fe_space%constraint_dofs_coefficients%get(pos)
          call this%recursive_matrix_assembly( test_fe, trial_fe, k, j, weight*a_ij, matrix, array)
        end do
      end if
  end if
end subroutine hp_adaptive_fe_face_iterator_recursive_matrix_assembly

recursive subroutine hp_adaptive_fe_face_iterator_recursive_vector_assembly(this, test_fe, i, v_i, array)
  implicit none
  class(hp_adaptive_fe_face_iterator_t), intent(in)    :: this
  type(hp_adaptive_fe_iterator_t)      , intent(in)    :: test_fe
  integer(ip)                          , intent(in)    :: i
  real(rp)                             , intent(in)    :: v_i
  type(serial_scalar_array_t)          , intent(inout) :: array
  
  integer(ip) :: k, pos, spos, epos
  real(rp) :: weight
 
  if ( .not. test_fe%is_fixed_dof(i) ) then    
      call array%add(i,v_i)
  else ! i is as fixed DoF
    if ( .not. test_fe%is_strong_dirichlet_dof(i) ) then
        ! Traverse DoFs on which i depends on
        spos = test_fe%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(i))
        epos = test_fe%hp_adaptive_fe_space%ptr_constraint_dofs%get(abs(i)+1)-1
        do pos=spos, epos 
          k       = test_fe%hp_adaptive_fe_space%constraint_dofs_dependencies%get(pos)
          weight  = test_fe%hp_adaptive_fe_space%constraint_dofs_coefficients%get(pos)
          call this%recursive_vector_assembly( test_fe, k, weight*v_i, array )
        end do
     end if
  end if
end subroutine hp_adaptive_fe_face_iterator_recursive_vector_assembly

subroutine serial_hp_adaptive_fe_space_create_fe_vef_iterator ( this, fe_vef )
  implicit none
  class(serial_hp_adaptive_fe_space_t), target, intent(in)    :: this
  type(fe_vef_iterator_t)                     , intent(inout) :: fe_vef
  type(p4est_vef_iterator_t) :: vef
  call fe_vef%free()
  call fe_vef%create(this,vef)
end subroutine serial_hp_adaptive_fe_space_create_fe_vef_iterator

subroutine serial_hp_adaptive_fe_space_create_fe_iterator ( this, fe )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(in)  :: this
  class(fe_iterator_t), allocatable, intent(inout) :: fe
  integer(ip) :: istat
  call this%free_fe_iterator(fe)
  allocate(hp_adaptive_fe_iterator_t :: fe, stat=istat); check(istat==0)
  call fe%create(this)
end subroutine serial_hp_adaptive_fe_space_create_fe_iterator

subroutine serial_hp_adaptive_fe_space_create_fe_face_iterator ( this, fe_face )
  implicit none
  class(serial_hp_adaptive_fe_space_t), target     , intent(in)    :: this
  class(fe_face_iterator_t)           , allocatable, intent(inout) :: fe_face
  type(p4est_vef_iterator_t) :: vef
  integer(ip)                :: istat
  call this%free_fe_face_iterator(fe_face)
  allocate(hp_adaptive_fe_face_iterator_t :: fe_face, stat=istat); check(istat==0)
  call fe_face%create(this,vef)
end subroutine serial_hp_adaptive_fe_space_create_fe_face_iterator

subroutine shpafs_create_same_reference_fes_on_all_cells ( this,          &
                                                           triangulation, &
                                                           reference_fes, &
                                                           conditions )
  implicit none
  class(serial_hp_adaptive_fe_space_t)        , intent(inout) :: this
  class(base_static_triangulation_t), target  , intent(in)    :: triangulation
  type(p_reference_fe_t)                      , intent(in)    :: reference_fes(:)
  class(conditions_t)       , target, optional, intent(in)    :: conditions

  type(std_vector_integer_ip_t), pointer :: vef_lids_of_fe_faces

  call this%free()

  call this%set_triangulation(triangulation) 
  select type(triangulation)
  class is (p4est_serial_triangulation_t)
    this%p4est_triangulation => triangulation
  class default
    assert(.false.)
  end select
 
  call this%set_number_fields(size(reference_fes))
  call this%allocate_and_fill_reference_fes(reference_fes)
  call this%allocate_ref_fe_id_per_fe()
  call this%fill_ref_fe_id_per_fe_same_on_all_cells()
  call this%check_cell_vs_fe_topology_consistency()
  call this%allocate_and_fill_fe_space_type_per_field()
  vef_lids_of_fe_faces => this%get_vef_lids_of_fe_faces()
  call vef_lids_of_fe_faces%resize(triangulation%get_num_vefs())
  call this%fill_vef_lids_of_fe_faces()
  call this%allocate_and_init_ptr_lst_dofs()
  
  if ( present(conditions) ) call this%set_conditions(conditions)
  call this%allocate_and_init_at_strong_dirichlet_bound()
  call this%allocate_and_init_has_fixed_dofs()
  call this%set_up_strong_dirichlet_bcs()
  
end subroutine shpafs_create_same_reference_fes_on_all_cells 

subroutine shpafs_create_different_between_cells( this,          &
                                                  triangulation, &
                                                  reference_fes, &
                                                  set_ids_to_reference_fes, &
                                                  conditions )
  implicit none
  class(serial_hp_adaptive_fe_space_t)        , intent(inout) :: this
  class(base_static_triangulation_t), target  , intent(in)    :: triangulation
  type(p_reference_fe_t)                      , intent(in)    :: reference_fes(:)
  integer(ip)                                 , intent(in)    :: set_ids_to_reference_fes(:,:)
  class(conditions_t)       , target, optional, intent(in)    :: conditions

  type(std_vector_integer_ip_t), pointer :: vef_lids_of_fe_faces

  call this%free()

  call this%set_triangulation(triangulation) 
  select type(triangulation)
  class is (p4est_serial_triangulation_t)
    this%p4est_triangulation => triangulation
  class default
    assert(.false.)
  end select
  
  call this%set_number_fields(size(set_ids_to_reference_fes,1))
  call this%allocate_and_fill_reference_fes(reference_fes)
  call this%allocate_ref_fe_id_per_fe()
  call this%fill_ref_fe_id_per_fe_different_between_cells(set_ids_to_reference_fes)
  call this%check_cell_vs_fe_topology_consistency()
  call this%allocate_and_fill_fe_space_type_per_field()
  vef_lids_of_fe_faces => this%get_vef_lids_of_fe_faces()
  call vef_lids_of_fe_faces%resize(triangulation%get_num_vefs())
  call this%fill_vef_lids_of_fe_faces()
  call this%allocate_and_init_ptr_lst_dofs()
  
  if ( present(conditions) ) call this%set_conditions(conditions)
  call this%allocate_and_init_at_strong_dirichlet_bound()
  call this%allocate_and_init_has_fixed_dofs()
  call this%set_up_strong_dirichlet_bcs()
  
end subroutine shpafs_create_different_between_cells

subroutine serial_hp_adaptive_fe_space_free(this)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout)    :: this
  call this%serial_fe_space_t%free()
  this%number_fixed_dofs = -1
  call this%free_ptr_constraint_dofs()
  call this%free_constraint_dofs_dependencies()
  call this%free_constraint_dofs_coefficients()
  nullify(this%p4est_triangulation)
end subroutine serial_hp_adaptive_fe_space_free

subroutine shpafs_free_ptr_constraint_dofs( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  call this%ptr_constraint_dofs%free()
end subroutine shpafs_free_ptr_constraint_dofs

subroutine shpafs_free_constraint_dofs_dependencies( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  call this%constraint_dofs_dependencies%free()  
end subroutine shpafs_free_constraint_dofs_dependencies

subroutine shpafs_free_constraint_dofs_coefficients( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: istat, block_id
  call this%constraint_dofs_coefficients%free()   
end subroutine shpafs_free_constraint_dofs_coefficients

function shpafs_get_number_fixed_dofs(this)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(in) :: this 
  integer(ip) :: shpafs_get_number_fixed_dofs
  shpafs_get_number_fixed_dofs = this%number_fixed_dofs
end function shpafs_get_number_fixed_dofs

subroutine shpafs_set_up_strong_dirichlet_bcs( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  integer(ip) :: i
  
  call serial_fe_space_set_up_strong_dirichlet_bcs(this)
  
  ! Re-size to 0 to force re-initialization during second resize (to the actual/correct size)
  call this%ptr_constraint_dofs%resize(0)
  this%number_fixed_dofs = this%get_number_strong_dirichlet_dofs()
  call this%ptr_constraint_dofs%resize(this%number_fixed_dofs+1,0)
  call this%constraint_dofs_dependencies%resize(this%number_fixed_dofs)
  call this%constraint_dofs_coefficients%resize(this%number_fixed_dofs)
  
  ! DoFs which are subject to strong Dirichlet BC's depend on an artificial
  ! DoF with identifier equal to zero. The corresponding value in constraint_dofs_coefficients
  ! is the value of the DoF subject to strong Dirichlet BC's
  call this%ptr_constraint_dofs%set(1,1)
  do i=2, this%ptr_constraint_dofs%size()
    call this%ptr_constraint_dofs%set(i, &
                                      this%ptr_constraint_dofs%get(i-1)+1 )
    
    call this%constraint_dofs_dependencies%set(i-1,0)
  end do
  
end subroutine shpafs_set_up_strong_dirichlet_bcs


subroutine shpafs_update_fixed_dof_values(this, free_dof_values, fixed_dof_values)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(in)    :: this 
  class(vector_t)                     , intent(in)    :: free_dof_values
  type(serial_scalar_array_t)         , intent(inout) :: fixed_dof_values
  real(rp) :: dof_value(1)
  real(rp) :: alpha 
  real(rp), pointer :: fixed_dof_values_entries(:)
  integer(ip) :: dof_lid, i, j, spos, epos
  
  fixed_dof_values_entries => fixed_dof_values%get_entries()
  do i=1, this%number_fixed_dofs
    spos=this%ptr_constraint_dofs%get(i)
    epos=this%ptr_constraint_dofs%get(i+1)-1
    alpha=0.0_rp
    do j=spos, epos
      dof_lid=this%constraint_dofs_dependencies%get(j)
      if ( dof_lid==0 ) then ! Dirichlet DoF
         alpha = this%constraint_dofs_coefficients%get(j)
         assert (spos == epos)
         exit
      else
         if (dof_lid > 0) then ! Depends on a free DoF
            call free_dof_values%extract_subvector( iblock=1, size_indices=1, indices=[dof_lid], values=dof_value )
         else ! Depends on a strong Dirichlet DoF
            call fixed_dof_values%extract_subvector( iblock=1, size_indices=1, indices=[-dof_lid], values=dof_value )
         end if
         alpha = alpha + this%constraint_dofs_coefficients%get(j) * dof_value(1)
      end if
    end do
    fixed_dof_values_entries(i) = alpha
  end do
end subroutine shpafs_update_fixed_dof_values

subroutine shpafs_interpolate_dirichlet_values (this, conditions, time, fields_to_interpolate)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout)  :: this
  class(conditions_t)     , intent(in)     :: conditions
  real(rp)       , optional   , intent(in) :: time
  integer(ip)    , optional   , intent(in) :: fields_to_interpolate(:)
  call serial_fe_space_interpolate_dirichlet_values(this, conditions, time, fields_to_interpolate)
  call this%transfer_dirichlet_to_constraint_dof_coefficients()
end subroutine shpafs_interpolate_dirichlet_values 

subroutine shpafs_fill_vef_lids_of_fe_faces ( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  class(base_static_triangulation_t), pointer     :: triangulation
  type(std_vector_integer_ip_t)     , pointer     :: vef_lids_of_fe_faces
  type(fe_vef_iterator_t)                         :: fe_vef
  integer(ip)                                     :: face_lid
  face_lid = 0
  triangulation => this%get_triangulation()
  vef_lids_of_fe_faces => this%get_vef_lids_of_fe_faces()
  call this%create_fe_vef_iterator(fe_vef)
  do while ( .not. fe_vef%has_finished() )
    if ( is_fe_face() ) then
      face_lid = face_lid + 1
      if ( face_lid <= vef_lids_of_fe_faces%size() ) then
        call fe_vef%set_vef_lid_of_fe_face(face_lid)
      else
        call fe_vef%push_back_vef_lid_of_fe_face()
      end if
    end if
    call fe_vef%next()
  end do
  call vef_lids_of_fe_faces%resize(face_lid)
  call this%free_fe_vef_iterator(fe_vef)
  
  contains
  
    function is_fe_face()
      implicit none
      logical :: is_fe_face
      is_fe_face = fe_vef%is_face() .and. ( ( .not. fe_vef%is_proper() ) .or.                         &
                                            ( fe_vef%is_proper() .and. ( fe_vef%is_at_boundary() .or. &
                                                                         fe_vef%get_num_cells_around() > 1 ) ) )
    end function is_fe_face
  
end subroutine shpafs_fill_vef_lids_of_fe_faces

subroutine serial_hp_adaptive_fe_space_fill_dof_info( this, block_layout )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  type(block_layout_t), target        , intent(inout) :: block_layout
  logical :: perform_numbering
  integer(ip) :: block_id, field_id
  type(block_layout_t), pointer :: p_block_layout
  
  p_block_layout => this%get_block_layout()
  perform_numbering = .not. associated(p_block_layout) 
  if (.not. perform_numbering) perform_numbering = .not. (p_block_layout == block_layout)
  
  if ( perform_numbering ) then
    call this%set_block_layout(block_layout)
  
    this%number_fixed_dofs = this%get_number_strong_dirichlet_dofs()
  
    ! Initialize number DoFs per field
    call this%allocate_number_dofs_per_field()
    do field_id=1, this%get_number_fields()
      call this%set_field_number_dofs(field_id, 0)
    end do
  
    ! Initialize number DoFs per block
    do block_id=1, this%get_number_blocks()
      call this%set_block_number_dofs(block_id, 0)
    end do
  
    ! Generate field-wise/block-wise global DoF identifiers
    do field_id = 1, this%get_number_fields()
      call this%fill_elem2dof_and_count_dofs( field_id )
    end do
  
    call this%setup_hanging_node_constraints()
  end if  
end subroutine serial_hp_adaptive_fe_space_fill_dof_info

subroutine serial_hp_adaptive_fe_space_fill_elem2dof_and_count_dofs( this, field_id ) 
  implicit none
  ! Parameters
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip)                         ,  intent(in)   :: field_id

  ! Local variables
  integer(ip) :: ivef, vef_lid
  integer(ip) :: iblock, init_dof_block, current_dof_block, previous_dof_block
  integer(ip) :: init_fixed_dof, current_fixed_dof, previous_fixed_dof
  integer(ip), allocatable :: visited_proper_vef_to_fe_map(:,:)
  integer(ip), allocatable :: visited_improper_vef_to_fe_map(:,:)
  
  class(fe_iterator_t) , allocatable :: fe, source_fe, coarser_fe
  type(fe_vef_iterator_t) :: vef
  integer(ip), pointer :: field_blocks(:)
  integer(ip), pointer :: fe_space_type_per_field(:)

  logical :: all_improper_cells_around_void, is_owner
  integer(ip) :: source_cell_id
  integer(ip) :: source_vef_lid  
  integer(ip) :: icell_improper_around

  field_blocks            => this%get_field_blocks()
  fe_space_type_per_field => this%get_fe_space_type()
  iblock            = field_blocks(field_id)
  init_dof_block    = this%get_block_number_dofs(iblock)
  current_dof_block = init_dof_block
  
  init_fixed_dof    = this%number_fixed_dofs
  current_fixed_dof = init_fixed_dof

  call this%create_fe_iterator(fe)
  if ( fe_space_type_per_field(field_id) == fe_space_type_cg ) then
     call memalloc ( 2, this%p4est_triangulation%get_num_proper_vefs(), visited_proper_vef_to_fe_map  ,  __FILE__, __LINE__ )
     call memalloc ( 2, this%p4est_triangulation%get_num_proper_vefs(), visited_improper_vef_to_fe_map,  __FILE__, __LINE__ )
     visited_proper_vef_to_fe_map = -1
     visited_improper_vef_to_fe_map = -1
     
     call this%create_fe_vef_iterator(vef)
     call this%create_fe_iterator(source_fe)
     call this%create_fe_iterator(coarser_fe)
     do while ( .not. fe%has_finished())
        if ( fe%is_local() ) then
           call fe%fill_own_dofs ( field_id, current_dof_block )
           do ivef = 1, fe%get_num_vefs()
              call fe%get_vef(ivef,vef)
              
              all_improper_cells_around_void=.true.
              do icell_improper_around=1, vef%get_num_improper_cells_around()
                call vef%get_improper_cell_around(icell_improper_around,coarser_fe)
                if (.not. coarser_fe%is_void(field_id)) then 
                   all_improper_cells_around_void=.false.
                   exit
                end if
              end do
              
              if ( vef%is_proper() .or. all_improper_cells_around_void ) then

                 vef_lid = abs(fe%get_vef_lid(ivef))
                 is_owner = .false.
                 if ( vef%is_proper()) then
                   is_owner = ( visited_proper_vef_to_fe_map   ( 1, vef_lid ) == -1 )
                 else
                   is_owner = ( visited_improper_vef_to_fe_map ( 1, vef_lid ) == -1 )
                 end if

                 if ( is_owner ) then
                    previous_dof_block = current_dof_block
                    call fe%fill_own_dofs_on_vef ( ivef, field_id, current_dof_block, free_dofs_loop=.true.  )
                    if (previous_dof_block < current_dof_block) then
                      if ( vef%is_proper()) then
                        visited_proper_vef_to_fe_map ( 1, vef_lid ) = fe%get_lid()
                        visited_proper_vef_to_fe_map ( 2, vef_lid ) = ivef
                      else
                        visited_improper_vef_to_fe_map ( 1, vef_lid ) = fe%get_lid()
                        visited_improper_vef_to_fe_map ( 2, vef_lid ) = ivef
                      end if
                    end if
                 else 
                    if ( vef%is_proper()) then
                      source_cell_id = visited_proper_vef_to_fe_map(1,vef_lid)
                      source_vef_lid = visited_proper_vef_to_fe_map(2,vef_lid)
                    else
                      source_cell_id = visited_improper_vef_to_fe_map(1,vef_lid)
                      source_vef_lid = visited_improper_vef_to_fe_map(2,vef_lid)
                    end if
                    call source_fe%set_lid(source_cell_id)
                    call fe%fill_own_dofs_on_vef_from_source_fe ( ivef, source_fe, source_vef_lid, field_id ) 
                 end if
              else 
                 assert ( fe%get_vef_lid(ivef) < 0 )
                 vef_lid = abs(fe%get_vef_lid(ivef))
                 if ( visited_improper_vef_to_fe_map ( 1, vef_lid ) == -1 ) then
                    previous_fixed_dof = current_fixed_dof
                    call fe%fill_own_dofs_on_vef ( ivef, field_id, current_fixed_dof, free_dofs_loop=.false.  )
                    if (previous_fixed_dof < current_fixed_dof) then
                      visited_improper_vef_to_fe_map ( 1, vef_lid ) = fe%get_lid()
                      visited_improper_vef_to_fe_map ( 2, vef_lid ) = ivef
                    end if
                 else 
                    call source_fe%set_lid(visited_improper_vef_to_fe_map(1,vef_lid))
                    call fe%fill_own_dofs_on_vef_from_source_fe ( ivef, &
                         source_fe, &
                         visited_improper_vef_to_fe_map(2,vef_lid), &
                         field_id) 
                 end if
              end if   
           end do
          end if 
          call fe%determine_has_fixed_dofs(field_id)
          call fe%next()
        end do
        call this%free_fe_iterator(source_fe)
        call this%free_fe_iterator(coarser_fe)
        call this%free_fe_vef_iterator(vef)
        call memfree ( visited_proper_vef_to_fe_map  ,  __FILE__, __LINE__ )
        call memfree ( visited_improper_vef_to_fe_map,  __FILE__, __LINE__ )
  else    
     ! TODO: this code is a verbatim copy of the one of its parent.
     !       we should better split the parent into additional TBPs
     !       to avoid code replication
     do while ( .not. fe%has_finished())
        if ( fe%is_local() ) then
           call fe%fill_own_dofs ( field_id, current_dof_block )
        end if
        call fe%next()
     end do
  end if
  call this%free_fe_iterator(fe)

  call this%set_field_number_dofs(field_id,current_dof_block - init_dof_block)
  call this%set_block_number_dofs(iblock, this%get_block_number_dofs(iblock) + & 
                                          this%get_field_number_dofs(field_id))
  this%number_fixed_dofs = this%number_fixed_dofs + current_fixed_dof - init_fixed_dof
end subroutine serial_hp_adaptive_fe_space_fill_elem2dof_and_count_dofs

subroutine shpafs_setup_hanging_node_constraints ( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  class(fe_iterator_t), allocatable :: fe
  class(fe_iterator_t), allocatable :: coarser_fe
  type(fe_vef_iterator_t) :: fe_vef, coarser_vef
  type(list_iterator_t) :: fe_own_dofs_on_vef_iterator
  type(list_iterator_t) :: fe_dofs_on_vef_iterator 
  type(list_iterator_t) :: coarser_fe_dofs_on_vef_iterator
  integer(ip) :: improper_vef_lid
  integer(ip) :: block_id, field_id
  class(reference_fe_t), pointer :: reference_fe, coarser_reference_fe
  type(i1p_t), allocatable :: elem2dof(:), coarser_fe_elem2dof(:)
  integer(ip) :: istat, i
  integer(ip) :: improper_dof_lid
  integer(ip) :: improper_vef_ivef, coarser_fe_ivef, coarse_fe_subvef
  integer(ip), pointer :: field_blocks(:)
  integer(ip) :: qpoint, ishape_fe, ishape_coarser_fe
  type(interpolation_t), pointer :: h_refinement_interpolation
  integer(ip), pointer :: h_refinement_subedge_permutation(:,:,:)
  integer(ip), pointer :: h_refinement_subface_permutation(:,:,:)
  real(rp) :: coefficient
  integer(ip) :: num_cell_vertices, num_cell_edges, num_cell_faces

  field_blocks => this%get_field_blocks()
  
  call this%ptr_constraint_dofs%resize(this%number_fixed_dofs+1,0)
  
  ! Transform header to length strong Dirichlet DoFs
  do improper_dof_lid=1, this%get_number_strong_dirichlet_dofs()
    call this%ptr_constraint_dofs%set(improper_dof_lid+1, 1 )
  end do

  call this%create_fe_iterator(fe)
  call this%create_fe_iterator(coarser_fe)
  call this%create_fe_vef_iterator(fe_vef)
  call this%create_fe_vef_iterator(coarser_vef)
  
  reference_fe => fe%get_reference_fe(1)
  num_cell_vertices = reference_fe%get_number_vertices()
  if  (this%p4est_triangulation%get_num_dimensions() == 3) then
    num_cell_edges = reference_fe%get_number_n_faces_of_dimension(1)
  else
    num_cell_edges = 0
  end if  
  num_cell_faces = reference_fe%get_number_faces()
  
  allocate(elem2dof(this%get_number_fields()), stat=istat); check(istat==0);
  allocate(coarser_fe_elem2dof(this%get_number_fields()), stat=istat); check(istat==0);

  ! Computation of constraints     
  do improper_vef_lid = 1, this%p4est_triangulation%get_num_improper_vefs()
     ! Retrieve all data related to the current improper vef 
     ! and one of the cells that owns it
     call fe_vef%set_lid(-improper_vef_lid)
     call fe_vef%get_cell_around(1,fe)
     improper_vef_ivef = fe%find_lpos_vef_lid(fe_vef%get_lid())
     call fe%get_elem2dof(elem2dof)

     ! Retrieve all data related to the first improper cell around current improper vef
     call fe_vef%get_improper_cell_around(1,coarser_fe)
     call coarser_fe%get_elem2dof(coarser_fe_elem2dof)
     coarser_fe_ivef = fe_vef%get_improper_cell_around_ivef()

     do field_id=1, this%get_number_fields()
        reference_fe => fe%get_reference_fe(field_id)
        coarser_reference_fe  => coarser_fe%get_reference_fe(field_id)
        block_id = field_blocks(field_id)
        fe_own_dofs_on_vef_iterator = reference_fe%create_own_dofs_on_n_face_iterator(improper_vef_ivef)
        do while (.not. fe_own_dofs_on_vef_iterator%is_upper_bound() )
           improper_dof_lid = elem2dof(field_id)%p(fe_own_dofs_on_vef_iterator%get_current())

           if ( .not. fe%is_fixed_dof(improper_dof_lid) ) exit

           improper_dof_lid = abs(improper_dof_lid)    

           coarser_fe_dofs_on_vef_iterator = coarser_reference_fe%create_dofs_on_n_face_iterator(coarser_fe_ivef)
           do while (.not. coarser_fe_dofs_on_vef_iterator %is_upper_bound() )
              call this%ptr_constraint_dofs%set(improper_dof_lid+1, &
                                                this%ptr_constraint_dofs%get(improper_dof_lid+1)+1 )
              call coarser_fe_dofs_on_vef_iterator%next()
           end do
           call fe_own_dofs_on_vef_iterator%next() 
        end do
     end do
  end do
  
  
  call this%p4est_triangulation%std_vector_transform_length_to_header(this%ptr_constraint_dofs)
  call this%constraint_dofs_dependencies%resize(this%ptr_constraint_dofs%get(this%ptr_constraint_dofs%size())-1)
  call this%constraint_dofs_coefficients%resize(this%ptr_constraint_dofs%get(this%ptr_constraint_dofs%size())-1)
  
  
  do improper_dof_lid=1, this%get_number_strong_dirichlet_dofs()
    call this%ptr_constraint_dofs%set(improper_dof_lid, this%ptr_constraint_dofs%get(improper_dof_lid+1) )
  end do
  
  
   ! Computation of constraints     
  do improper_vef_lid = 1, this%p4est_triangulation%get_num_improper_vefs()
     ! Retrieve all data related to the current improper vef 
     ! and one of the cells that owns it
     call fe_vef%set_lid(-improper_vef_lid)
     call fe_vef%get_cell_around(1,fe)
     improper_vef_ivef = fe%find_lpos_vef_lid(fe_vef%get_lid())
     call fe%get_elem2dof(elem2dof)

     mcheck( this%p4est_triangulation%get_num_dimensions()==2 , 'The following code only valid for 2d cases' )
     ! Here we need to find the first non-void coarser_fe in this field. In 2d, the first non-void coarser_fe is always the first
     ! one (if it exists). But this is not true in 3d. Thus, taking the first coarser_fe in the next line is correct only in 2d.

     ! Retrieve all data related to the first improper cell around current improper vef
     call fe_vef%get_improper_cell_around(1,coarser_fe)
     call coarser_fe%get_elem2dof(coarser_fe_elem2dof)
     coarser_fe_ivef  = fe_vef%get_improper_cell_around_ivef()
     coarse_fe_subvef = fe_vef%get_improper_cell_around_subvef()
     call coarser_fe%get_vef(coarser_fe_ivef,coarser_vef)

     do field_id=1, this%get_number_fields()
        if ( coarser_fe%is_void(field_id)) cycle

        reference_fe => fe%get_reference_fe(field_id)
        coarser_reference_fe => coarser_fe%get_reference_fe(field_id)
        
        select type(coarser_reference_fe)
        type is (hex_lagrangian_reference_fe_t)
           h_refinement_subedge_permutation => coarser_reference_fe%get_h_refinement_subedge_permutation()
           h_refinement_subface_permutation => coarser_reference_fe%get_h_refinement_subface_permutation()
        class default
          assert(.false.)
        end select
        
        block_id = field_blocks(field_id)
        fe_own_dofs_on_vef_iterator = reference_fe%create_own_dofs_on_n_face_iterator(improper_vef_ivef)
        fe_dofs_on_vef_iterator = reference_fe%create_dofs_on_n_face_iterator(improper_vef_ivef)
        do while (.not. fe_own_dofs_on_vef_iterator%is_upper_bound() )
           ishape_fe = fe_own_dofs_on_vef_iterator%get_current()
           improper_dof_lid = elem2dof(field_id)%p(ishape_fe)
           assert ( fe%is_fixed_dof(improper_dof_lid) )
           improper_dof_lid = abs(improper_dof_lid)   
           
           call fe_dofs_on_vef_iterator%begin() 
           do while (.not. fe_dofs_on_vef_iterator%is_upper_bound() )
             if ( fe_dofs_on_vef_iterator%get_current() == ishape_fe ) exit
             call fe_dofs_on_vef_iterator%next() 
           end do
           assert (.not. fe_dofs_on_vef_iterator%is_upper_bound() )
           
           if ( fe_vef%get_dimension() == 0 ) then ! vef is a corner (2D/3D)
              if ( coarser_vef%get_dimension()  == 1 .and. this%p4est_triangulation%get_num_dimensions() == 3) then
                 !qpoint = h_refinement_subedge_permutation(coarser_fe_ivef,num_subedges_per_edge,1)
              else
                 ! TODO: 2 is equal to num_subfaces_per_face ... only working in 2D!!!
                 qpoint = h_refinement_subface_permutation(coarser_fe_ivef-num_cell_vertices,2,1)
              end if
           else if ( fe_vef%get_dimension()  == 1 .and. this%p4est_triangulation%get_num_dimensions() == 3 ) then ! vef is an edge (only 3D)
              qpoint = h_refinement_subedge_permutation(coarser_fe_ivef-num_cell_vertices, &
                                                        coarse_fe_subvef, &
                                                        fe_dofs_on_vef_iterator%get_distance_to_lower_bound())
           else if (fe_vef%get_dimension() == this%p4est_triangulation%get_num_dimensions()-1) then ! vef is a face (2D/3D)
              qpoint = h_refinement_subface_permutation(coarser_fe_ivef-num_cell_vertices-num_cell_edges, &
                                                        coarse_fe_subvef, &
                                                        fe_dofs_on_vef_iterator%get_distance_to_lower_bound())
           end if
           
           coarser_fe_dofs_on_vef_iterator = coarser_reference_fe%create_dofs_on_n_face_iterator(coarser_fe_ivef)
           do while (.not. coarser_fe_dofs_on_vef_iterator %is_upper_bound() )
              
              ishape_coarser_fe = coarser_fe_dofs_on_vef_iterator%get_current() 
              call this%constraint_dofs_dependencies%set(this%ptr_constraint_dofs%get(improper_dof_lid),&
                                                         coarser_fe_elem2dof(field_id)%p(ishape_coarser_fe))

              ! Evaluate coefficient and push_back into the corresponding std_vector data structure
              select type(coarser_reference_fe)
                type is (hex_lagrangian_reference_fe_t)
                call coarser_reference_fe%get_h_refinement_coefficient(ishape_fe,ishape_coarser_fe,qpoint,coefficient) 
              class default
                assert(.false.)
              end select
              call this%constraint_dofs_coefficients%set(this%ptr_constraint_dofs%get(improper_dof_lid),&
                                                         coefficient)

              call this%ptr_constraint_dofs%set(improper_dof_lid, &
                                                this%ptr_constraint_dofs%get(improper_dof_lid)+1 )
              call coarser_fe_dofs_on_vef_iterator%next()
           end do
           call fe_own_dofs_on_vef_iterator%next() 
        end do
     end do
  end do
  
  do i=this%ptr_constraint_dofs%size(),2,-1
    call this%ptr_constraint_dofs%set(i, this%ptr_constraint_dofs%get(i-1))
  end do
  call this%ptr_constraint_dofs%set(1,1)
  
  ! TODO: shrink_to_fit this%constraint_dofs_dependencies and this%constraint_dofs_coefficients
  call this%free_fe_iterator(coarser_fe)
  call this%free_fe_iterator(fe)
  call this%free_fe_vef_iterator(fe_vef)
  call this%free_fe_vef_iterator(coarser_vef)
end subroutine shpafs_setup_hanging_node_constraints 

subroutine shpafs_transfer_dirichlet_to_constraint_dof_coefficients(this)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: improper_dof_lid
  type(serial_scalar_array_t), pointer :: strong_dirichlet_values
  real(rp), pointer :: strong_dirichlet_values_entries(:)
  strong_dirichlet_values         => this%get_strong_dirichlet_values()
  strong_dirichlet_values_entries => strong_dirichlet_values%get_entries()
  do improper_dof_lid=1, this%get_number_strong_dirichlet_dofs()
     call this%constraint_dofs_coefficients%set(improper_dof_lid, strong_dirichlet_values_entries(improper_dof_lid) )
  end do
end subroutine shpafs_transfer_dirichlet_to_constraint_dof_coefficients

subroutine shpafs_transfer_dirichlet_to_fe_space(this,fixed_dof_values)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  type(serial_scalar_array_t)         , intent(in)    :: fixed_dof_values
  type(serial_scalar_array_t), pointer     :: strong_dirichlet_values
  real(rp)                   , pointer     :: strong_dirichlet_values_entries(:)
  integer(ip)                              :: i, num_strong_dirichlet_dofs
  integer(ip)                , allocatable :: indices(:)
  num_strong_dirichlet_dofs = this%get_number_strong_dirichlet_dofs()
  call memalloc(num_strong_dirichlet_dofs,indices,__FILE__,__LINE__)
  indices = (/ (i, i=1,num_strong_dirichlet_dofs) /)
  strong_dirichlet_values => this%get_strong_dirichlet_values()
  strong_dirichlet_values_entries => strong_dirichlet_values%get_entries()
  call fixed_dof_values%extract_subvector( 1, num_strong_dirichlet_dofs, &
                                           indices, strong_dirichlet_values_entries )
  call this%transfer_dirichlet_to_constraint_dof_coefficients()
  call memfree(indices,__FILE__,__LINE__)
end subroutine shpafs_transfer_dirichlet_to_fe_space

subroutine shpafs_project_ref_fe_id_per_fe(this)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  class(base_static_triangulation_t), pointer     :: triangulation
  type(std_vector_integer_ip_t)     , pointer     :: p4est_refinement_and_coarsening_flags
  class(fe_iterator_t)              , allocatable :: new_fe
  class(reference_fe_t)             , pointer     :: reference_fe
  integer(ip)                       , allocatable :: old_ref_fe_id_per_fe(:,:)
  integer(ip)                                     :: old_num_cells
  integer(ip)                                     :: num_children_per_cell
  integer(ip)                                     :: subcell_id, old_cell_lid, new_cell_lid
  integer(ip)                                     :: current_old_cell_lid, current_new_cell_lid
  integer(ip)                                     :: field_id
  integer(ip)                                     :: transformation_flag
  integer(ip)                                     :: old_reference_fe_id
  
  triangulation => this%get_triangulation()
  select type(triangulation)
  class is (p4est_serial_triangulation_t)
    p4est_refinement_and_coarsening_flags => triangulation%get_p4est_refinement_and_coarsening_flags()
  class default
    assert(.false.)
  end select
  
  call this%move_alloc_ref_fe_id_per_fe_out(old_ref_fe_id_per_fe)
  call this%allocate_ref_fe_id_per_fe()
  
  old_num_cells = size(old_ref_fe_id_per_fe,2)
  
  call this%create_fe_iterator(new_fe)
  
  reference_fe => new_fe%get_reference_fe_geo()
  num_children_per_cell = reference_fe%get_number_n_faces_of_dimension(0)
  
  old_cell_lid = 1
  new_cell_lid = 1
  do while ( old_cell_lid .le. old_num_cells )
    transformation_flag = p4est_refinement_and_coarsening_flags%get(old_cell_lid)
    call new_fe%set_lid(new_cell_lid)
    do field_id = 1,this%get_number_fields()
      current_old_cell_lid = old_cell_lid
      current_new_cell_lid = new_cell_lid
      old_reference_fe_id  = old_ref_fe_id_per_fe(field_id,current_old_cell_lid)
      if ( transformation_flag == do_nothing ) then
        call new_fe%set_reference_fe_id(field_id,old_reference_fe_id)
        current_new_cell_lid = current_new_cell_lid + 1
      else if ( transformation_flag == refinement ) then
        do subcell_id = 0,num_children_per_cell-1
          call new_fe%set_reference_fe_id(field_id,old_reference_fe_id)
          current_new_cell_lid = current_new_cell_lid + 1
          call new_fe%set_lid(current_new_cell_lid)
        end do
      else if ( transformation_flag == coarsening ) then
        do subcell_id = 1,num_children_per_cell-1
          current_old_cell_lid = current_old_cell_lid + 1
          if ( old_reference_fe_id /= old_ref_fe_id_per_fe(field_id,current_old_cell_lid) ) then
            massert(.false.,'Coarsened subcells do not have the same reference FE id')
          end if
        end do
        call new_fe%set_reference_fe_id(field_id,old_reference_fe_id)
        current_new_cell_lid = current_new_cell_lid + 1
      else
        massert(.false.,'Unrecognised refinement and coarsening flag')
      end if
    end do
    old_cell_lid = current_old_cell_lid
    new_cell_lid = current_new_cell_lid
    old_cell_lid = old_cell_lid + 1
  end do
  
  massert ( new_cell_lid - 1 == this%p4est_triangulation%get_num_cells(), 'Loop in old cells failed to visit all new cells' )
  
  call this%free_fe_iterator(new_fe)
  call memfree(old_ref_fe_id_per_fe,__FILE__,__LINE__)
  
end subroutine shpafs_project_ref_fe_id_per_fe

subroutine shpafs_project_fe_integration_arrays(this)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  ! if allocated, FE integration arrays not projected
end subroutine shpafs_project_fe_integration_arrays

subroutine shpafs_project_fe_face_integration_arrays(this)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  ! if allocated, FE face integration arrays not projected
end subroutine shpafs_project_fe_face_integration_arrays

subroutine serial_hp_adaptive_fe_space_refine_and_coarsen( this, fe_function )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  type(fe_function_t)                 , intent(inout) :: fe_function
  class(base_static_triangulation_t), pointer     :: triangulation
  type(fe_function_t)                             :: transformed_fe_function
  type(std_vector_integer_ip_t)     , pointer     :: p4est_refinement_and_coarsening_flags
  integer(ip)                       , allocatable :: old_ptr_dofs_per_fe(:,:)
  integer(ip)              , target , allocatable :: old_lst_dofs_lids(:)
  integer(ip)                       , pointer     :: old_field_elem2dof(:)
  integer(ip)                                     :: num_children_per_cell
  integer(ip)                                     :: transformation_flag
  integer(ip)                                     :: subcell_id, old_cell_lid, new_cell_lid
  integer(ip)                                     :: current_old_cell_lid, current_new_cell_lid
  integer(ip)                                     :: old_num_cells
  integer(ip)                                     :: field_id
  class(fe_iterator_t)              , allocatable :: new_fe
  real(rp)                          , allocatable :: old_nodal_values(:,:)
  real(rp)                          , allocatable :: new_nodal_values(:)
  class(reference_fe_t)             , pointer     :: reference_fe
  integer(ip)                                     :: number_nodes_field
  type(block_layout_t)              , pointer     :: block_layout
  
  triangulation => this%get_triangulation()
  select type(triangulation)
  class is (p4est_serial_triangulation_t)
    p4est_refinement_and_coarsening_flags => triangulation%get_p4est_refinement_and_coarsening_flags()
  class default
    assert(.false.)
  end select
  
  old_num_cells = p4est_refinement_and_coarsening_flags%size()
  
  call this%move_alloc_ptr_dofs_per_fe_out(old_ptr_dofs_per_fe)
  call this%move_alloc_lst_dofs_lids_out(old_lst_dofs_lids)
  massert ( old_num_cells == (size(old_ptr_dofs_per_fe,2)-1), 'Incorrect size of p4est_refinement_and_coarsening_flags' )
  
  call this%project_ref_fe_id_per_fe()
  !call this%check_cell_vs_fe_topology_consistency()
  call this%fill_vef_lids_of_fe_faces()
  call this%project_fe_integration_arrays()
  call this%project_fe_face_integration_arrays()
  call this%allocate_and_init_ptr_lst_dofs()
  call this%allocate_and_init_at_strong_dirichlet_bound()
  call this%allocate_and_init_has_fixed_dofs()
  call this%set_up_strong_dirichlet_bcs()
  
  ! Force that a new DoF numbering is generated for the refined/coarsened triangulation
  block_layout => this%get_block_layout()
  call this%nullify_block_layout()
  call this%fill_dof_info(block_layout)
  
  call this%create_fe_iterator(new_fe)
  reference_fe => new_fe%get_reference_fe_geo()
  num_children_per_cell = reference_fe%get_number_n_faces_of_dimension(0)
  call memalloc(num_children_per_cell, &
                this%get_max_number_shape_functions(),old_nodal_values,__FILE__,__LINE__)
  call memalloc(this%get_max_number_shape_functions(),new_nodal_values,__FILE__,__LINE__)
  
  call transformed_fe_function%create(this)
  
  old_cell_lid = 1
  new_cell_lid = 1
  do while ( old_cell_lid .le. old_num_cells )
    transformation_flag = p4est_refinement_and_coarsening_flags%get(old_cell_lid)
    call new_fe%set_lid(new_cell_lid)
    do field_id = 1,this%get_number_fields()
      current_old_cell_lid = old_cell_lid
      current_new_cell_lid = new_cell_lid
      reference_fe => new_fe%get_reference_fe(field_id) ! Only h-adaptivity
      number_nodes_field = reference_fe%get_number_shape_functions()
      old_field_elem2dof => get_field_elem2dof()
      call fe_function%gather_nodal_values( field_id,                & 
                                            old_field_elem2dof,      &
                                            number_nodes_field,      & 
                                            this%get_field_blocks(), &
                                            old_nodal_values(1,1:number_nodes_field) )
      if ( transformation_flag == do_nothing ) then
        call transformed_fe_function%insert_nodal_values( new_fe,   &
                                                          field_id, &
                                                          old_nodal_values(1,1:number_nodes_field) )
        current_new_cell_lid = current_new_cell_lid + 1
      else if ( transformation_flag == refinement ) then
        do subcell_id = 0,num_children_per_cell-1
          select type(reference_fe)
          type is (hex_lagrangian_reference_fe_t)
            call reference_fe%interpolate_nodal_values_on_subcell( subcell_id,                               & 
                                                                   old_nodal_values(1,1:number_nodes_field), &
                                                                   new_nodal_values(1:number_nodes_field) )
            call transformed_fe_function%insert_nodal_values( new_fe,   &
                                                              field_id, &
                                                              new_nodal_values(1:number_nodes_field) )
          type is (void_reference_fe_t)
            ! Do nothing
          class default
            assert(.false.)
          end select
          current_new_cell_lid = current_new_cell_lid + 1
          call new_fe%set_lid(current_new_cell_lid)
        end do
      else if ( transformation_flag == coarsening ) then
        do subcell_id = 1,num_children_per_cell-1
          current_old_cell_lid = current_old_cell_lid + 1
          old_field_elem2dof => get_field_elem2dof()
          call fe_function%gather_nodal_values( field_id,                & 
                                                old_field_elem2dof,      &
                                                number_nodes_field,      & 
                                                this%get_field_blocks(), &
                                                old_nodal_values(subcell_id+1,1:number_nodes_field) )
        end do
        select type(reference_fe)
        type is (hex_lagrangian_reference_fe_t)
          call reference_fe%project_nodal_values_on_cell( old_nodal_values(:,1:number_nodes_field), &
                                                          new_nodal_values(1:number_nodes_field) )
          call transformed_fe_function%insert_nodal_values( new_fe,   &
                                                            field_id, &
                                                            new_nodal_values(1:number_nodes_field) )
        type is (void_reference_fe_t)
          ! Do nothing
        class default
          assert(.false.)
        end select
        current_new_cell_lid = current_new_cell_lid + 1
      else
        massert(.false.,'Unrecognised refinement and coarsening flag')
      end if
    end do
    old_cell_lid = current_old_cell_lid
    new_cell_lid = current_new_cell_lid
    old_cell_lid = old_cell_lid + 1
  end do

  massert ( new_cell_lid - 1 == this%p4est_triangulation%get_num_cells(), 'Loop in old cells failed to visit all new cells' )
  
  call fe_function%create(this)
  fe_function = transformed_fe_function
  
  call this%transfer_dirichlet_to_fe_space( fe_function%get_fixed_dof_values() )
  
  select type(triangulation)
  class is (p4est_serial_triangulation_t)
    call triangulation%clear_refinement_and_coarsening_flags()
  class default
    assert(.false.)
  end select
  
  call this%free_fe_iterator(new_fe)
  call transformed_fe_function%free()
  call memfree(old_nodal_values,__FILE__,__LINE__)
  call memfree(new_nodal_values,__FILE__,__LINE__)
  call memfree(old_ptr_dofs_per_fe,__FILE__,__LINE__)
  call memfree(old_lst_dofs_lids,__FILE__,__LINE__)
  
contains
  
  function get_field_elem2dof()
    implicit none
    integer(ip), pointer     :: get_field_elem2dof(:)
    integer(ip)              :: spos, epos
    spos = old_ptr_dofs_per_fe(field_id,current_old_cell_lid)
    if ( field_id == this%get_number_fields() ) then
      epos = old_ptr_dofs_per_fe(1,current_old_cell_lid+1)-1
    else
      epos = old_ptr_dofs_per_fe(field_id+1,current_old_cell_lid)-1
    end if
    get_field_elem2dof => old_lst_dofs_lids(spos:epos)    
  end function get_field_elem2dof
  
end subroutine serial_hp_adaptive_fe_space_refine_and_coarsen


end module hp_adaptive_fe_space_names
