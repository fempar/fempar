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
  use environment_names
  use conditions_names
  use std_vector_integer_ip_names
  use std_vector_real_rp_names
  use list_types_names
  
  
  implicit none
# include "debug.i90"
  private

  ! TODO: to move member variables of serial_hp_adaptive_fe_space_t to serial_fe_space_t
  type, extends(serial_fe_space_t) :: serial_hp_adaptive_fe_space_t
     private
     !list_t :: constraints 
     ! Hanging and strong Dirichlet data per block
     ! ptr_constraint_dofs         : pointer to constraints
     ! l1                          : constraint DOFs dependencies (0 for independent term)
     ! constraint_dofs_coefficients: constraint DoFs coefficients (also independent term)
     ! u_fixed = sum u_dep w_dep + c
     type(std_vector_integer_ip_t), allocatable :: ptr_constraint_dofs(:)
     type(std_vector_integer_ip_t), allocatable :: constraint_dofs_dependencies(:)
     type(std_vector_real_rp_t)   , allocatable :: constraint_dofs_coefficients(:)
     
     type(p4est_serial_triangulation_t), pointer :: p4est_triangulation =>  NULL()
     
     integer(ip), allocatable :: number_fixed_dofs_per_field(:)
     integer(ip), allocatable :: number_fixed_dofs_per_block(:)
     
   contains
     procedure                            :: create_fe_vef_iterator                                 => serial_hp_adaptive_fe_space_create_fe_vef_iterator
     
     
     procedure                            :: serial_fe_space_create_same_reference_fes_on_all_cells => shpafs_create_same_reference_fes_on_all_cells 
     procedure                            :: free                                                   => serial_hp_adaptive_fe_space_free
     
     
     procedure                            :: fill_dof_info                                          => serial_hp_adaptive_fe_space_fill_dof_info
     procedure                  , private :: fill_elem2dof_and_count_dofs                           => serial_hp_adaptive_fe_space_fill_elem2dof_and_count_dofs
     procedure                            :: allocate_number_fixed_dofs_per_block                   => shpafs_allocate_number_fixed_dofs_per_block
     procedure                            :: free_number_fixed_dofs_per_block                       => shpafs_free_number_fixed_dofs_per_block
     procedure                            :: allocate_number_fixed_dofs_per_field                   => shpafs_allocate_number_fixed_dofs_per_field
     procedure                            :: free_number_fixed_dofs_per_field                       => shpafs_free_number_fixed_dofs_per_field
     
     procedure                            :: setup_hanging_node_constraints                         => shpafs_setup_hanging_node_constraints
     procedure                            :: allocate_ptr_constraint_dofs                           => shpafs_allocate_ptr_constraint_dofs
     procedure                            :: free_ptr_constraint_dofs                               => shpafs_free_ptr_constraint_dofs
     procedure                            :: allocate_constraint_dofs_dependencies                  => shpafs_allocate_constraint_dofs_dependencies
     procedure                            :: free_constraint_dofs_dependencies                      => shpafs_free_constraint_dofs_dependencies
     procedure                            :: allocate_constraint_dofs_coefficients                  => shpafs_allocate_constraint_dofs_coefficients
     procedure                            :: free_constraint_dofs_coefficients                      => shpafs_free_constraint_dofs_coefficients
     
     ! Getters & Setters
     procedure                            :: get_field_number_fixed_dofs                            => serial_hp_adaptive_fe_space_get_field_number_fixed_dofs
     procedure                            :: set_field_number_fixed_dofs                            => serial_hp_adaptive_fe_space_set_field_number_fixed_dofs
     procedure                            :: get_block_number_fixed_dofs                            => serial_hp_adaptive_fe_space_get_block_number_fixed_dofs
     procedure                            :: set_block_number_fixed_dofs                            => serial_hp_adaptive_fe_space_set_block_number_fixed_dofs
 end type serial_hp_adaptive_fe_space_t  
 
 public :: serial_hp_adaptive_fe_space_t
 
contains

subroutine serial_hp_adaptive_fe_space_create_fe_vef_iterator ( this, fe_vef )
  implicit none
  class(serial_hp_adaptive_fe_space_t), target, intent(in)    :: this
  type(fe_vef_iterator_t)                     , intent(inout) :: fe_vef
  type(p4est_vef_iterator_t) :: vef
  call fe_vef%free()
  call fe_vef%create(this,vef)
end subroutine serial_hp_adaptive_fe_space_create_fe_vef_iterator

subroutine shpafs_create_same_reference_fes_on_all_cells ( this, &
                                                           triangulation, &
                                                           conditions, &
                                                           reference_fes, &
                                                           field_blocks, &
                                                           field_coupling )
  implicit none
  class(serial_hp_adaptive_fe_space_t)        , intent(inout) :: this
  class(base_static_triangulation_t), target  , intent(in)    :: triangulation
  class(conditions_t)                         , intent(in)    :: conditions
  type(p_reference_fe_t)                      ,  intent(in)   :: reference_fes(:)
  integer(ip)                       , optional, intent(in)    :: field_blocks(:)
  logical                           , optional, intent(in)    :: field_coupling(:,:)

  integer(ip) :: i, istat, jfield, ifield

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
  call this%allocate_and_fill_field_blocks_and_coupling(field_blocks, field_coupling)
  call this%allocate_ref_fe_id_per_fe()
  call this%fill_ref_fe_id_per_fe_same_on_all_cells()
  call this%check_cell_vs_fe_topology_consistency()
  call this%allocate_and_fill_fe_space_type_per_field()
  call this%allocate_and_init_ptr_lst_dofs()
  call this%allocate_and_init_at_strong_dirichlet_bound()
  call this%set_up_strong_dirichlet_bcs( conditions )
end subroutine shpafs_create_same_reference_fes_on_all_cells 

subroutine serial_hp_adaptive_fe_space_free(this)
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout)    :: this
  call this%serial_fe_space_t%free()
  call this%free_number_fixed_dofs_per_block()
  call this%free_number_fixed_dofs_per_field()
  call this%free_ptr_constraint_dofs()
  call this%free_constraint_dofs_dependencies()
  call this%free_constraint_dofs_coefficients()
  nullify(this%p4est_triangulation)
end subroutine   

subroutine shpafs_allocate_number_fixed_dofs_per_block( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  if ( .not. allocated(this%number_fixed_dofs_per_block) ) & 
     call memalloc( this%get_number_blocks(), this%number_fixed_dofs_per_block, __FILE__, __LINE__ )
  assert ( size(this%number_fixed_dofs_per_block) == this%get_number_blocks() ) 
end subroutine shpafs_allocate_number_fixed_dofs_per_block

subroutine shpafs_free_number_fixed_dofs_per_block( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  if ( allocated(this%number_fixed_dofs_per_block) ) & 
     call memfree( this%number_fixed_dofs_per_block, __FILE__, __LINE__ )
end subroutine shpafs_free_number_fixed_dofs_per_block

subroutine shpafs_allocate_number_fixed_dofs_per_field( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  if ( .not. allocated(this%number_fixed_dofs_per_field) ) & 
     call memalloc( this%get_number_fields(), this%number_fixed_dofs_per_field, __FILE__, __LINE__ )
  assert ( size(this%number_fixed_dofs_per_field) == this%get_number_fields() )
end subroutine shpafs_allocate_number_fixed_dofs_per_field

subroutine shpafs_free_number_fixed_dofs_per_field( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  if ( allocated(this%number_fixed_dofs_per_field) ) & 
     call memfree( this%number_fixed_dofs_per_field, __FILE__, __LINE__ )
end subroutine shpafs_free_number_fixed_dofs_per_field

subroutine shpafs_allocate_ptr_constraint_dofs( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: istat
  
  if ( .not. allocated(this%ptr_constraint_dofs) ) then
    allocate(this%ptr_constraint_dofs(this%get_number_blocks()), stat=istat ); assert (istat == 0)
  end if
  assert ( size (this%ptr_constraint_dofs) == this%get_number_blocks() )
end subroutine shpafs_allocate_ptr_constraint_dofs

subroutine shpafs_free_ptr_constraint_dofs( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: istat, block_id
  
  if ( allocated(this%ptr_constraint_dofs) ) then
    do block_id=1, this%get_number_blocks()
      call this%ptr_constraint_dofs(block_id)%free()
    end do  
    deallocate(this%ptr_constraint_dofs, stat=istat ); assert (istat == 0)
  end if
end subroutine shpafs_free_ptr_constraint_dofs

subroutine shpafs_allocate_constraint_dofs_dependencies( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: istat
  
  if ( .not. allocated(this%constraint_dofs_dependencies) ) then
    allocate(this%constraint_dofs_dependencies(this%get_number_blocks()), stat=istat ); assert (istat == 0)
  end if
  assert ( size (this%constraint_dofs_dependencies) == this%get_number_blocks() )
end subroutine shpafs_allocate_constraint_dofs_dependencies

subroutine shpafs_free_constraint_dofs_dependencies( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: istat, block_id
  
  if ( allocated(this%constraint_dofs_dependencies) ) then
    do block_id=1, this%get_number_blocks()
      call this%constraint_dofs_dependencies(block_id)%free()
    end do  
    deallocate(this%constraint_dofs_dependencies, stat=istat ); assert (istat == 0)
  end if
end subroutine shpafs_free_constraint_dofs_dependencies

subroutine shpafs_allocate_constraint_dofs_coefficients( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: istat
  
  if ( .not. allocated(this%constraint_dofs_coefficients) ) then
    allocate(this%constraint_dofs_coefficients(this%get_number_blocks()), stat=istat ); assert (istat == 0)
  end if
  assert ( size (this%constraint_dofs_coefficients) == this%get_number_blocks() )
end subroutine shpafs_allocate_constraint_dofs_coefficients

subroutine shpafs_free_constraint_dofs_coefficients( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: istat, block_id
  
  if ( allocated(this%constraint_dofs_coefficients) ) then
    do block_id=1, this%get_number_blocks()
      call this%constraint_dofs_coefficients(block_id)%free()
    end do  
    deallocate(this%constraint_dofs_coefficients, stat=istat ); assert (istat == 0)
  end if
end subroutine shpafs_free_constraint_dofs_coefficients

subroutine serial_hp_adaptive_fe_space_fill_dof_info( this )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip) :: block_id, field_id

  call this%allocate_number_dofs_per_block()
  call this%allocate_number_fixed_dofs_per_block()
  do block_id=1, this%get_number_blocks()
    call this%set_block_number_dofs(block_id, 0)
    call this%set_block_number_fixed_dofs(block_id, 0)
  end do
  
  call this%allocate_number_dofs_per_field()
  call this%allocate_number_fixed_dofs_per_field()
  do field_id=1, this%get_number_fields()
    call this%set_field_number_dofs(field_id, 0)
    call this%set_field_number_fixed_dofs(field_id, 0)
  end do
  
  do field_id = 1, this%get_number_fields()
     call this%fill_elem2dof_and_count_dofs( field_id )
  end do
  
  call this%setup_hanging_node_constraints()
end subroutine serial_hp_adaptive_fe_space_fill_dof_info

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
  integer(ip) :: qpoint, ishape
  
  call this%allocate_ptr_constraint_dofs()
  call this%allocate_constraint_dofs_dependencies()
  call this%allocate_constraint_dofs_coefficients()

  field_blocks => this%get_field_blocks()
  do block_id = 1, this%get_number_blocks()
     ! Re-size to 0 to force re-initialization during second resize (to the actual/correct size)
     call this%ptr_constraint_dofs(block_id)%resize(0)
     call this%ptr_constraint_dofs(block_id)%resize(this%get_block_number_fixed_dofs(block_id)+1,0)
  end do

  call this%create_fe_iterator(fe)
  call this%create_fe_iterator(coarser_fe)
  call this%create_fe_vef_iterator(fe_vef)
  call this%create_fe_vef_iterator(coarser_vef)

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
           assert ( improper_dof_lid < 0 )
           improper_dof_lid = abs(improper_dof_lid)    

           coarser_fe_dofs_on_vef_iterator = coarser_reference_fe%create_dofs_on_n_face_iterator(coarser_fe_ivef)
           do while (.not. coarser_fe_dofs_on_vef_iterator %is_upper_bound() )
              call this%ptr_constraint_dofs(block_id)%set(improper_dof_lid+1, &
                                                          this%ptr_constraint_dofs(block_id)%get(improper_dof_lid+1)+1 )
              call coarser_fe_dofs_on_vef_iterator%next()
           end do
           call fe_own_dofs_on_vef_iterator%next() 
        end do
     end do
  end do
  
  do block_id=1, this%get_number_blocks()
     call this%p4est_triangulation%std_vector_transform_length_to_header(this%ptr_constraint_dofs(block_id))
     call this%constraint_dofs_dependencies(block_id)%resize(this%ptr_constraint_dofs(block_id)%get(this%ptr_constraint_dofs(block_id)%size())-1)
     call this%constraint_dofs_coefficients(block_id)%resize(this%ptr_constraint_dofs(block_id)%get(this%ptr_constraint_dofs(block_id)%size())-1)
  end do
  
  
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
     coarser_fe_ivef  = fe_vef%get_improper_cell_around_ivef()
     coarse_fe_subvef = fe_vef%get_improper_cell_around_subvef()
     call coarser_fe%get_vef(coarser_fe_ivef,coarser_vef)

     do field_id=1, this%get_number_fields()
        reference_fe => fe%get_reference_fe(field_id)
        coarser_reference_fe => coarser_fe%get_reference_fe(field_id)
        block_id = field_blocks(field_id)
        fe_own_dofs_on_vef_iterator = reference_fe%create_own_dofs_on_n_face_iterator(improper_vef_ivef)
        fe_dofs_on_vef_iterator = reference_fe%create_dofs_on_n_face_iterator(improper_vef_ivef)
        do while (.not. fe_own_dofs_on_vef_iterator%is_upper_bound() )
           improper_dof_lid = elem2dof(field_id)%p(fe_own_dofs_on_vef_iterator%get_current())
           assert ( improper_dof_lid < 0 )
           improper_dof_lid = abs(improper_dof_lid)   
           
           call fe_dofs_on_vef_iterator%begin() 
           do while (.not. fe_dofs_on_vef_iterator%is_upper_bound() )
             if ( fe_dofs_on_vef_iterator%get_current() == fe_own_dofs_on_vef_iterator%get_current() ) exit
             call fe_dofs_on_vef_iterator%next() 
           end do
           assert (.not. fe_dofs_on_vef_iterator%is_upper_bound() )
           
           if ( fe_vef%get_dimension() == 0 ) then ! vef is a corner (2D/3D)
              if ( coarser_vef%get_dimension()  == 1 .and. this%p4est_triangulation%get_num_dimensions() == 3) then
                 !qpoint = coarser_fe%h_refinement_edge_permutation(coarse_fe_ivef,num_subedges_per_edge,1)
              else
                 !qpoint = coarser_fe%h_refinement_face_permutation(coarse_fe_ivef,num_subface_per_face,1)
              end if
           else if ( fe_vef%get_dimension()  == 1 .and. this%p4est_triangulation%get_num_dimensions() == 3 ) then ! vef is an edge (only 3D)
              !qpoint = coarser_fe%h_refinement_edge_permutation(coarse_fe_ivef,coarse_fe_subvef,fe_dofs_on_vef_iterator%get_distance_to_lower_bound())
           else if (fe_vef%get_dimension() == this%p4est_triangulation%get_num_dimensions()-1) then ! vef is a face (2D/3D)
              ! Extract those quadrature points identifiers which lay on ivef+subface
              !qpoint = coarser_fe%h_refinement_face_permutation(coarse_fe_ivef,coarse_fe_subvef,fe_dofs_on_vef_iterator%get_distance_to_lower_bound())
           end if
            
           coarser_fe_dofs_on_vef_iterator = coarser_reference_fe%create_dofs_on_n_face_iterator(coarser_fe_ivef)
           do while (.not. coarser_fe_dofs_on_vef_iterator %is_upper_bound() )
              
              ishape = coarser_fe_dofs_on_vef_iterator%get_current() 
              call this%constraint_dofs_dependencies(block_id)%set(this%ptr_constraint_dofs(block_id)%get(improper_dof_lid),&
                                                                   coarser_fe_elem2dof(field_id)%p(ishape))
              
              ! Evaluate coefficient and push_back into the corresponding std_vector data structure
              !call coarse_reference_fe%get_value(ishape, qpoint, interpolation_h_refinement, coefficient) 
              call this%constraint_dofs_coefficients(block_id)%set(this%ptr_constraint_dofs(block_id)%get(improper_dof_lid),&
                                                                    0.0_rp)

              call this%ptr_constraint_dofs(block_id)%set(improper_dof_lid, &
                                                          this%ptr_constraint_dofs(block_id)%get(improper_dof_lid)+1 )
              call coarser_fe_dofs_on_vef_iterator%next()
           end do
           call fe_own_dofs_on_vef_iterator%next() 
        end do
     end do
  end do
  
  do block_id=1, this%get_number_blocks()
     do i=this%ptr_constraint_dofs(block_id)%size(),2,-1
       call this%ptr_constraint_dofs(block_id)%set(i, this%ptr_constraint_dofs(block_id)%get(i-1))
     end do
     call this%ptr_constraint_dofs(block_id)%set(1,1) 
  end do
  
  ! TODO: shrink_to_fit this%constraint_dofs_dependencies and this%constraint_dofs_coefficients
  call this%free_fe_iterator(coarser_fe)
  call this%free_fe_iterator(fe)
  call this%free_fe_vef_iterator(fe_vef)
  call this%free_fe_vef_iterator(coarser_vef)
end subroutine shpafs_setup_hanging_node_constraints 

subroutine serial_hp_adaptive_fe_space_fill_elem2dof_and_count_dofs( this, field_id ) 
  implicit none
  ! Parameters
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this 
  integer(ip)                         ,  intent(in)   :: field_id

  ! Local variables
  integer(ip) :: ivef, vef_lid, ielem
  integer(ip) :: iblock, init_dof_block, current_dof_block, previous_dof_block
  integer(ip) :: init_fixed_dof_block, current_fixed_dof_block, previous_fixed_dof_block
  integer(ip), allocatable :: visited_proper_vef_to_fe_map(:,:)
  integer(ip), allocatable :: visited_improper_vef_to_fe_map(:,:)
  
  class(fe_iterator_t) , allocatable :: fe, source_fe
  type(fe_vef_iterator_t) :: vef
  integer(ip), pointer :: field_blocks(:)
  integer(ip), pointer :: fe_space_type_per_field(:)
  
  field_blocks            => this%get_field_blocks()
  fe_space_type_per_field => this%get_fe_space_type()
  iblock            = field_blocks(field_id)
  init_dof_block    = this%get_block_number_dofs(iblock)
  current_dof_block = init_dof_block
  
  init_fixed_dof_block    = this%get_block_number_fixed_dofs(iblock)
  current_fixed_dof_block = init_fixed_dof_block

  call this%create_fe_iterator(fe)
  if ( fe_space_type_per_field(field_id) == fe_space_type_cg ) then
     call memalloc ( 2, this%p4est_triangulation%get_num_proper_vefs(), visited_proper_vef_to_fe_map  ,  __FILE__, __LINE__ )
     call memalloc ( 2, this%p4est_triangulation%get_num_proper_vefs(), visited_improper_vef_to_fe_map,  __FILE__, __LINE__ )
     visited_proper_vef_to_fe_map = -1
     visited_improper_vef_to_fe_map = -1
     
     call this%create_fe_vef_iterator(vef)
     call this%create_fe_iterator(source_fe)
     do while ( .not. fe%has_finished())
        if ( fe%is_local() ) then
           call fe%fill_own_dofs ( field_id, current_dof_block )
           do ivef = 1, fe%get_num_vefs()
              call fe%get_vef(ivef,vef)
              if ( vef%is_proper() ) then
                 vef_lid = fe%get_vef_lid(ivef)
                 assert ( vef_lid > 0 )
                 if ( visited_proper_vef_to_fe_map ( 1, vef_lid ) == -1 ) then
                    previous_dof_block = current_dof_block
                    call fe%fill_own_dofs_on_vef ( ivef, field_id, current_dof_block, free_dofs_loop=.true.  )
                    if (previous_dof_block < current_dof_block) then
                      visited_proper_vef_to_fe_map ( 1, vef_lid ) = fe%get_lid()
                      visited_proper_vef_to_fe_map ( 2, vef_lid ) = ivef
                    end if
                 else 
                    call source_fe%set_lid(visited_proper_vef_to_fe_map(1,vef_lid))
                    call fe%fill_own_dofs_on_vef_from_source_fe ( ivef, &
                         source_fe, &
                         visited_proper_vef_to_fe_map(2,vef_lid), &
                         field_id) 
                 end if
              else 
                 assert ( fe%get_vef_lid(ivef) < 0 )
                 vef_lid = abs(fe%get_vef_lid(ivef))
                 if ( visited_improper_vef_to_fe_map ( 1, vef_lid ) == -1 ) then
                    previous_fixed_dof_block = current_fixed_dof_block
                    call fe%fill_own_dofs_on_vef ( ivef, field_id, current_fixed_dof_block, free_dofs_loop=.false.  )
                    if (previous_fixed_dof_block < current_fixed_dof_block) then
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
          call fe%next()
        end do
        call this%free_fe_iterator(source_fe)
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
  call this%set_field_number_fixed_dofs(field_id,current_fixed_dof_block - init_fixed_dof_block)
  call this%set_block_number_dofs(iblock, this%get_block_number_dofs(iblock) + & 
                                          this%get_field_number_dofs(field_id))
  call this%set_block_number_fixed_dofs(iblock,this%get_block_number_fixed_dofs(iblock) + & 
                                               this%get_field_number_fixed_dofs(field_id))
end subroutine serial_hp_adaptive_fe_space_fill_elem2dof_and_count_dofs


! Returns the number of fixed DoFs associated to field with identifier field_id
function serial_hp_adaptive_fe_space_get_field_number_fixed_dofs( this, field_id )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(in) :: this
  integer(ip)          , intent(in) :: field_id
  integer(ip)                       :: serial_hp_adaptive_fe_space_get_field_number_fixed_dofs
  class(environment_t), pointer  :: environment
  environment => this%get_environment()
  serial_hp_adaptive_fe_space_get_field_number_fixed_dofs = 0
  if ( environment%am_i_l1_task() ) then
     assert ( field_id >=1 .and. field_id <= this%get_number_fields() ) 
     serial_hp_adaptive_fe_space_get_field_number_fixed_dofs = this%number_fixed_dofs_per_field(field_id)
  end if
end function serial_hp_adaptive_fe_space_get_field_number_fixed_dofs

subroutine serial_hp_adaptive_fe_space_set_field_number_fixed_dofs( this, field_id, field_number_fixed_dofs )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  integer(ip)           , intent(in)    :: field_id
  integer(ip)           , intent(in)    :: field_number_fixed_dofs
  class(environment_t), pointer  :: environment
  environment => this%get_environment()
  if ( environment%am_i_l1_task() ) then
     assert ( field_id >=1 .and. field_id <= this%get_number_fields() )
     this%number_fixed_dofs_per_field(field_id) = field_number_fixed_dofs
  end if
end subroutine serial_hp_adaptive_fe_space_set_field_number_fixed_dofs

! Returns the number of fixed DoFs associated to block with identifier block_id
function serial_hp_adaptive_fe_space_get_block_number_fixed_dofs ( this, block_id )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(in) :: this
  integer(ip)          , intent(in) :: block_id
  integer(ip)                       :: serial_hp_adaptive_fe_space_get_block_number_fixed_dofs
  class(environment_t), pointer  :: environment
  environment => this%get_environment()
  serial_hp_adaptive_fe_space_get_block_number_fixed_dofs = 0
  if ( environment%am_i_l1_task() ) then
     assert ( block_id >=1 .and. block_id <= this%get_number_blocks() ) 
     serial_hp_adaptive_fe_space_get_block_number_fixed_dofs  = this%number_fixed_dofs_per_block(block_id)
  end if
end function serial_hp_adaptive_fe_space_get_block_number_fixed_dofs

subroutine serial_hp_adaptive_fe_space_set_block_number_fixed_dofs ( this, block_id, block_number_fixed_dofs )
  implicit none
  class(serial_hp_adaptive_fe_space_t), intent(inout) :: this
  integer(ip)           , intent(in)    :: block_id
  integer(ip)                           :: block_number_fixed_dofs
  class(environment_t), pointer         :: environment
  environment => this%get_environment()
  if ( environment%am_i_l1_task() ) then
     assert ( block_id >=1 .and. block_id <= this%get_number_blocks() )
     this%number_fixed_dofs_per_block(block_id) = block_number_fixed_dofs
  end if
end subroutine serial_hp_adaptive_fe_space_set_block_number_fixed_dofs



!! Assembly of local matrices for hp-adaptivity
!subroutine fe_iterator_assemble(this,elmat,elvec,matrix_array_assembler)
!  implicit none
!  class(fe_iterator_t)            , intent(in)   :: this
!  real(rp)                        , intent(in)   :: elmat(:,:)
!  real(rp)                        , intent(in)   :: elvec(:)
!  class(matrix_array_assembler_t), intent(inout) :: matrix_array_assembler
!  
!  call this%get_elem2dof(elem2dof)
!  
!  matrix => matrix_array_assembler%get_matrix()
!  array  => matrix_array_assembler%get_array()
!  
!  ! TODO: if only proper vefs on current cell call 
!  !      call matrix_array_assembler%assembly( number_fields, 
!  !                                            num_dofs_per_field, 
!  !                                            elem2dof, 
!  !                                            field_blocks, 
!  !                                            field_coupling, 
!  !                                            elmat, 
!  !                                            elvec )

! ielmat=0
! do ife_space=1, number_fe_spaces
!   iblock = field_blocks(ife_space)
!   call this%get_elem2dof(ife_space,ielem2dof)
!   jelmat=0
!   do jfe_space=1, number_fe_spaces
!     if ((field_coupling(ife_space,jfe_space))) then
!        jblock = field_blocks(jfe_space)
!        call this%get_elem2dof(jfe_space,jelem2dof)
!        
!        scalar_matrix => matrix
!        select type(matrix)
!          type is (block_sparse_matrix_t)
!            scalar_matrix => matrix%get_block(iblock,jblock)
!          type is (par_block_sparse_matrix_t)
!            scalar_matrix => matrix%get_block(iblock,jblock)
!        end select
!        
!        ! TODO: repeat select type block with scalar_array 
!        do i=1, this%get_num_dofs(ife_space)
!         do j=1, this%get_num_dofs(jfe_space)
!           call this%recursive_assembly( ielemdof%(i), jelem2dof(j), elmat(ielmat+i,jelmat+j), elvec(jelmat+j), scalar_matrix, scalar_array)
!         end do
!        end do      
!     end if
!     jelmat=jelmat+this%get_num_dofs(jfe_space)
!   end do
!   ielmat=ielmat+this%get_num_dofs(ife_space)
! end do
!end subroutine   
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





end module hp_adaptive_fe_space_names
