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
  ! Serial modules
  use types_names
  use p4est_serial_triangulation_names
  use reference_fe_names
  
  
  implicit none
# include "debug.i90"
  private

!  ! TODO: to move this stuff to serial_fe_space_t
!  type :: serial_hp_adaptive_fe_space_t  
!     private     
!     ! Reference FE container
!     integer(ip)                                 :: reference_fes_size
!     type(p_reference_fe_t)        , allocatable :: reference_fes(:)
!    
!     ! DoF identifiers associated to each FE and field within FE
!     integer(ip)                   , allocatable :: ptr_dofs_per_fe(:,:) ! (number_fields, number_fes+1)
!     integer(ip)                   , allocatable :: lst_dofs_lids(:)
!     
!     
!     list_t :: constraints 
!     ! Hanging and strong Dirichlet
!     ! p  : pointer to constraints
!     ! l1 : constraint dependency DOFs (0 for independent term)
!     ! l2 : constraint coefficients (also independent term)
!     ! u_fixed = sum u_dep w_dep + c
!     
!     integer(ip)                                 :: number_strong_dirichlet_dofs
!     type(serial_scalar_array_t)                 :: strong_dirichlet_values
!     logical                       , allocatable :: at_strong_dirichlet_boundary_per_fe(:,:)
!     
!     type(p4est_serial_triangulation_t), pointer :: triangulation =>  NULL()
!   contains
!     procedure                           :: fill_dof_info                                => serial_fe_space_fill_dof_info
!     procedure                 , private :: fill_elem2dof_and_count_dofs                 => serial_fe_space_fill_elem2dof_and_count_dofs
!     procedure                 , private :: renumber_dofs_block                          => serial_fe_space_renumber_dofs_block
! end type serial_hp_adaptive_fe_space_t  
! 
! 
!contains

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
!        k       = this%constraint_dependency_dofs%get(pos)
!        if (k>0) then ! If current DoF on which i depends is subject to Dirichlet BC's
!          weight  = this%constraint_coefficient%get(pos)
!          call this%recursive_assembly( k, j, weight*a_ij, v_j, matrix_array_assembler )
!        end if  
!      end do
!  else
!    if ( j < 0 ) then
!      ! Traverse DoFs on which j depends on
!      do pos=this%ptr_constraint(abs(j)), this%ptr_constraint(abs(j)+1)-1
!        k       = this%constraint_dependency_dofs%get(pos)
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









!end subroutine

!subroutine fill_dof_info

!  class(hp_adaptive_fe_iterator_t), allocatable :: fe
!  !   type(std_vector_integer_ip_t) :: proper_dofs_on_vefs
!  !   type(std_vector_integer_ip_t) :: improper_dofs_on_vefs
!  !   integer(ip), allocatable :: ptr_proper_dofs_on_vefs(:)
!  !   integer(ip), allocatable :: ptr_improper_dofs_on_vefs(:)
!  integer(ip), allocatable :: cell_owner_proper_vefs(:,:)
!  integer(ip), allocatable :: cell_owner_improper_vefs(:,:)



!  iblock            = this%field_blocks(field_id)
!  init_dof_block    = this%number_dofs_per_block(iblock)
!  current_dof_block = init_dof_block


!  call memalloc ( 2, this%triangulation%get_num_vefs(), visited_vef_to_fe_map,  __FILE__, __LINE__ )
!  visited_vef_to_fe_map = -1

!  call this%create_fe_iterator(source_fe)
!  call this%create_fe_iterator(fe)
!  do while ( .not. fe%has_finished())
!     if ( fe%is_local() ) then
!        call fe%fill_own_dofs ( field_id, current_dof_block )
!        do ivef = 1, fe%get_num_vefs()
!           call fe%get_vef(ivef,vef)
!           if ( vef%is_proper() ) then
!              vef_lid = fe%get_vef_lid(ivef)
!              if ( visited_vef_to_fe_map ( 1, vef_lid ) == -1 ) then
!                 visited_vef_to_fe_map ( 1, vef_lid ) = fe%get_lid()
!                 visited_vef_to_fe_map ( 2, vef_lid ) = ivef
!                 call fe%fill_own_dofs_on_vef ( ivef, field_id, current_dof_block  ) 
!              else 
!                 call source_fe%set_lid(visited_vef_to_fe_map(1,vef_lid))
!                 call fe%fill_own_dofs_on_vef_from_source_fe ( ivef, &
!                      source_fe, &
!                      visited_vef_to_fe_map(2,vef_lid), &
!                      field_id) 
!              end if
!           else    
!              vef_lid = abs(fe%get_vef_lid(ivef))
!              if ( visited_improper_vef_to_fe_map ( 1, vef_lid ) == -1 ) then
!                 visited_improper_vef_to_fe_map ( 1, vef_lid ) = fe%get_lid()
!                 visited_improper_vef_to_fe_map ( 2, vef_lid ) = ivef
!                 call fe%fill_own_dofs_on_vef ( ivef, field_id, improper_current_dof_block, .false.  ) 
!              else 
!                 call source_fe%set_lid(visited_improper_vef_to_fe_map(1,vef_lid))
!                 call fe%fill_own_dofs_on_vef_from_source_fe ( ivef, &
!                      source_fe, &
!                      visited_improper_vef_to_fe_map(2,vef_lid), &
!                      field_id) 
!              end if
!           end do
!        end if
!        call fe%next()
!     end do

!end subroutine fill_dof_info

!subroutine setup_hanging_nodes_constraints(visited_improper_vef_to_fe_map)
!    
!     ! TODO: REPEAT this loop  for each block!

!     ! Computation of constraints     
!     do improper_vef_lid = 1, this%triangulation%get_num_improper_vefs()
!        improper_vef_owner_cell_lid = visited_improper_vef_to_fe_map( 1, improper_vef_lid )
!        improper_vef_ivef = visited_improper_vef_to_fe_map( 2, improper_vef_lid )
!        call fe%set_lid(improper_vef_owner_cell_lid)
!        call fe%get_vef(improper_vef_ivef,vef)

!        ! Extract improper DoFs on current improper VEF????
!        call fe%get_elem2dof(elem2dof)
!        reference_fe => fe%get_reference_fe(field_id)
!        fe_own_dofs_on_vef_iterator = reference_fe%create_own_dofs_on_n_face_iterator(improper_vef_ivef)
!        fe_dofs_on_vef_iterator = reference_fe%create_dofs_on_n_face_iterator(improper_vef_ivef)
!        do while (.not. fe_own_dofs_on_vef_iterator%is_upper_bound() )
!           improper_dof_lid = elem2dof(fe_own_dofs_on_vef_iterator%current())
!           assert ( improper_dof_lid < 0 )
!           improper_dof_lid  = abs(improper_dof_lid)
!           
!           call fe_dofs_on_vef_iterator%init() 
!           do while (.not. fe_dofs_on_vef_iterator%is_upper_bound() )
!             if ( fe_dofs_on_vef_iterator%current() == fe_own_dofs_on_vef_iterator%current() ) exit
!             call fe_dofs_on_vef_iterator%next() 
!           end do
!           assert (.not. fe_dofs_on_vef_iterator%is_upper_bound() )

!           ! Extract all DOFs in the (closed) VEF of the coarser fe
!           call vef%get_improper_cell_around( 1, coarser_fe )
!           coarser_fe_ivef             = vef%get_ivef_improper_cell_around() 
!           coarse_fe_ivef_subentity_id = vef%get_ivef_subentity_id() 
!           call coarser_fe%get_vef(coarser_fe_ivef,coarser_vef)
!           if ( vef%get_dimension() == 0 ) then ! vef is a corner (2D/3D)
!              if ( coarse_vef%get_dimension()  == 1 .and. this%triangulation%get_num_dimensions() == 3) then
!                 qpoint = coarser_fe%h_refinement_edge_permutation(coarse_fe_ivef,num_subedges_per_edge,1)
!              else
!                 qpoint = coarser_fe%h_refinement_face_permutation(coarse_fe_ivef,num_subface_per_face,1)
!              end if
!           else if ( vef%get_dimension()  == 1 .and. this%triangulation%get_num_dimensions() == 3 ) then ! vef is an edge (only 3D)
!              qpoint = coarser_fe%h_refinement_edge_permutation(coarse_fe_ivef,coarse_fe_ivef_subentity_id,fe_dofs_on_vef_iterator%get_distance_to_lower_bound())
!           else (vef%get_dimension() == this%triangulation%get_num_dimensions()-1) then ! vef is a face (2D/3D)
!              coarser_fe_ivef_subface = vef%get_ivef_subface_improper_cell_around() 
!              ! Extract those quadrature points identifiers which lay on ivef+subface
!              qpoint = coarser_fe%h_refinement_face_permutation(coarse_fe_ivef,coarse_fe_ivef_subentity_id,fe_dofs_on_vef_iterator%get_distance_to_lower_bound())
!           end if
!           
!           call coarser_fe%get_elem2dof(coarser_fe_elem2dof)
!           coarser_reference_fe       => coarser_fe%get_reference_fe(field_id)
!           interpolation_h_refinement => coarser_fe%get_interpolation_h_refinement()
!           coarser_fe_dofs_on_vef_iterator = coarser_fe_reference_fe%create_dofs_on_n_face_iterator(coarser_fe_ivef)
!           do while (.not. coarser_fe_dofs_on_vef_iterator %is_upper_bound() )
!              this%ptr_constraint_dofs(improper_dof_lid+1) = this%ptr_constraint_dofs(improper_dof_lid+1) + 1
!              ishape = coarser_fe_dofs_on_vef_iterator%current() 
!              call this%constraint_dependency_dofs%push_back(coarser_fe_elem2dof(ishape))
!              ! Evaluate coefficient and push_back into the corresponding std_vector data structure
!              call coarse_reference_fe%get_value(ishape, qpoint, interpolation_h_refinement, coefficient) 
!              call this%constraint_coefficients%push_back(coefficient)
!              call coarser_fe_dofs_on_vef_iterator%next()
!           end do
!           call fe_own_dofs_on_vef_iterator%next()
!        end do
!     end do

!     ! Count to header transformation
!     this%ptr_constraint_dofs(1) = 1 
!     do improper_dof_lid = 1, this%get_num_improper_dof_lids(block_id)
!        this%ptr_constraint_dofs(improper_dof_lid+1) =  this%ptr_constraint_dofs(improper_dof_lid+1) + &
!             this%ptr_constraint_dofs(improper_dof_lid)
!     end do

!     ! TODO: shrink_to_fit this%constraint_dependency_dofs and this%constraint_coefficients
!     call this%free_fe_iterator(source_fe)
!     call this%free_fe_iterator(fe)
!end subroutine setup_hanging_nodes_constraints





end module hp_adaptive_fe_space_names
