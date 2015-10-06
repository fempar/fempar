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
module par_fe_space_names
  ! Serial Modules
  use types_names
  use memor_names
  use serial_fe_space_names
  use fe_space_types_names
  use problem_names
  use dof_descriptor_names
  use hash_table_names
  use finite_element_names

  ! Par Modules
  use par_environment_names
  use par_triangulation_names
  use par_element_exchange_names
  use par_conditions_names
  use blocks_dof_distribution_names
  
    ! Abstract modules
  use fe_space_names
  use assembler_names
  use matrix_array_assembler_names
  use matrix_names
  use array_names
  
  ! Concrete implementations 
  use par_scalar_matrix_array_assembler_names
  use par_scalar_matrix_names
  use par_block_matrix_names
  use par_scalar_array_names
  use par_block_array_names

#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private

  type, extends(fe_space_t) :: par_fe_space_t
     type(serial_fe_space_t) :: serial_fe_space
     type(par_triangulation_t), pointer :: p_trian => NULL()
	 
	 ! Extra data members that p_fe_space adds to fe_space (local portion)
	 type(blocks_dof_distribution_t) :: blocks_dof_distribution
	 integer(ip)                     :: num_interface_faces
     type(fe_face_t), allocatable    :: interface_faces(:)
  contains
     procedure :: create => par_fe_space_create
	 procedure :: free => par_fe_space_free
	 procedure :: print => par_fe_space_print
     procedure :: set_analytical_code => par_fe_space_set_analytical_code
	 procedure :: get_blocks_dof_distribution => par_fe_space_get_blocks_dof_distribution
	 
	 procedure, private :: par_fe_space_create_make_par_scalar_coefficient_matrix
     procedure, private :: par_fe_space_create_make_par_block_coefficient_matrix
	 generic :: make_coefficient_matrix => par_fe_space_create_make_par_scalar_coefficient_matrix, &
	                                       par_fe_space_create_make_par_block_coefficient_matrix
										   
	 procedure :: create_matrix_array_assembler         => par_fe_space_create_matrix_array_assembler
	 procedure :: symbolic_setup_matrix_array_assembler => par_fe_space_symbolic_setup_matrix_array_assembler
	 procedure :: volume_integral                       => par_fe_space_volume_integral 									   
  end type par_fe_space_t

  ! Types
  public :: par_fe_space_t

contains

 function par_fe_space_create_matrix_array_assembler(this,& 
											         diagonal_blocks_symmetric_storage,&
											         diagonal_blocks_symmetric,&
											         diagonal_blocks_sign)
    implicit none
	class(par_fe_space_t)          , intent(in) :: this
	logical                        , intent(in) :: diagonal_blocks_symmetric_storage(:)
    logical                        , intent(in) :: diagonal_blocks_symmetric(:)
	integer(ip)                    , intent(in) :: diagonal_blocks_sign(:)
	class(matrix_array_assembler_t), pointer    :: par_fe_space_create_matrix_array_assembler
	
	class(matrix_t), pointer :: matrix
	class(array_t) , pointer :: array
	
	assert ( size(diagonal_blocks_symmetric_storage) == this%dof_descriptor%nblocks )
	assert ( size(diagonal_blocks_symmetric) == this%dof_descriptor%nblocks )
	assert ( size(diagonal_blocks_sign) == this%dof_descriptor%nblocks )
	
	! 1. Select dynamically the type of class(matrix_array_assembler_t), class(matrix_t) and class(vector_t)
	! 2. Create class(matrix_t) and class(vector_t) accordingly to their dynamic type
	if (this%dof_descriptor%nblocks == 1) then
	  allocate ( par_scalar_matrix_array_assembler_t :: par_fe_space_create_matrix_array_assembler )
	  allocate ( par_scalar_matrix_t :: matrix )
	  allocate ( par_scalar_array_t  :: array )
	  select type(matrix)
        class is(par_scalar_matrix_t)
	      call matrix%create(diagonal_blocks_symmetric_storage(1),& 
							 diagonal_blocks_symmetric(1),& 
							 diagonal_blocks_sign(1), &
							 this%blocks_dof_distribution%get_block(1), &
						     this%p_trian%p_env)
	    class default
          check(.false.)
        end select 
	  select type(array)
        class is(par_scalar_array_t)
	      call array%create(this%blocks_dof_distribution%blocks(1), this%p_trian%p_env)
		  array%state = part_summed
	    class default
         check(.false.)
      end select 
	else
	  ! allocate ( par_block_matrix_array_assembler_t :: par_fe_space_create_matrix_array_assembler )
	  allocate ( par_block_matrix_t :: matrix )
	  allocate ( par_block_array_t  :: array )
	  check(.false.)
	  ! This itinerary still to be implemented ...
	end if
	call par_fe_space_create_matrix_array_assembler%set_matrix(matrix)
	call par_fe_space_create_matrix_array_assembler%set_array(array)
  end function par_fe_space_create_matrix_array_assembler
  
  subroutine par_fe_space_symbolic_setup_matrix_array_assembler(this,matrix_array_assembler)
	implicit none
	class(par_fe_space_t)           , intent(in)    :: this
	class(matrix_array_assembler_t) , intent(inout) :: matrix_array_assembler
	
    ! Polymorphic matrix 
	class(matrix_t), pointer :: matrix
	
	matrix => matrix_array_assembler%get_matrix()
    select type(matrix)
      class is(par_scalar_matrix_t)
	     if ( matrix%p_env%am_i_fine_task() ) then
            call setup_dof_graph_from_block_row_col_identifiers ( 1, 1, this%serial_fe_space, matrix%serial_scalar_matrix%graph )
		 end if	
	  class is(par_block_matrix_t)
	  class default
         check(.false.)
    end select
	
  end subroutine par_fe_space_symbolic_setup_matrix_array_assembler

  subroutine par_fe_space_volume_integral(this,approximations,assembler)
	implicit none
	class(par_fe_space_t)          , intent(in)    :: this
	class(p_discrete_integration_t), intent(in)    :: approximations(:) 
	class(assembler_t)             , intent(inout) :: assembler
	if ( this%p_trian%p_env%am_i_fine_task() ) then
	  call this%serial_fe_space%volume_integral(approximations,assembler)
	end if
  end subroutine par_fe_space_volume_integral

  subroutine par_fe_space_create_make_par_scalar_coefficient_matrix(this,symmetric_storage,is_symmetric,sign,par_scalar_matrix)
    implicit none
	class(par_fe_space_t)        , intent(in)  :: this
    logical                      , intent(in)  :: symmetric_storage
    logical                      , intent(in)  :: is_symmetric
	integer(ip)                  , intent(in)  :: sign
    type(par_scalar_matrix_t)    , intent(out) :: par_scalar_matrix
	
	assert ( this%dof_descriptor%nblocks == 1 ) 
	call par_scalar_matrix%create(symmetric_storage, & 
								  is_symmetric, &
								  sign, &
								  this%blocks_dof_distribution%get_block(1), &
								  this%p_trian%p_env)
	
	if ( this%p_trian%p_env%am_i_fine_task() ) then
	  call setup_dof_graph_from_block_row_col_identifiers ( 1, 1, this%serial_fe_space, par_scalar_matrix%serial_scalar_matrix%graph )
	end if  
	
	call par_scalar_matrix%allocate()
  end subroutine par_fe_space_create_make_par_scalar_coefficient_matrix
  
  subroutine par_fe_space_create_make_par_block_coefficient_matrix(this, & 
												  diagonal_blocks_symmetric_storage, & 
												  diagonal_blocks_symmetric, &
												  diagonal_blocks_sign,& 
												  par_block_matrix)
    implicit none
    class(par_fe_space_t)       , intent(in)  :: this
	logical                     , intent(in)  :: diagonal_blocks_symmetric_storage(this%dof_descriptor%nblocks)
    logical                     , intent(in)  :: diagonal_blocks_symmetric(this%dof_descriptor%nblocks)
    integer(ip)                 , intent(in)  :: diagonal_blocks_sign(this%dof_descriptor%nblocks)
	type(par_block_matrix_t)    , intent(out) :: par_block_matrix
	
	
  end subroutine par_fe_space_create_make_par_block_coefficient_matrix

  !*********************************************************************************
  ! This subroutine is intended to be called from a parallel driver, having as input
  ! a set of arrays of size number of local elements times number of global variables
  ! for continuity and order, and number of local elements for material and problem,
  ! together with some optional flags. The output of this subroutine is a fe_space
  ! with the required info on ghost elements.
  !*********************************************************************************
  subroutine par_fe_space_create ( this, p_trian, dof_descriptor, problem, p_cond, continuity, enable_face_integration,  &
                                    order, material,  time_steps_to_store, hierarchical_basis, &
                                    static_condensation, num_continuity )
    implicit none
    ! Dummy arguments
	class(par_fe_space_t)            , intent(inout) :: this  
    type(par_triangulation_t), target, intent(in)    :: p_trian
    type(dof_descriptor_t)           , intent(in)    :: dof_descriptor
    integer(ip)                      , intent(in)    :: problem(:)
    type(par_conditions_t)           , intent(in)    :: p_cond
    integer(ip)                      , intent(in)    :: continuity(:,:)
	logical                          , intent(in)    :: enable_face_integration(:,:)
    integer(ip)                      , intent(in)    :: order(:,:)
    integer(ip)                      , intent(in)    :: material(:)
    integer(ip)          , optional  , intent(in)    :: time_steps_to_store
    logical          , optional, intent(in)          :: hierarchical_basis
    logical          , optional, intent(in)          :: static_condensation
    integer(ip)          , optional, intent(in)      :: num_continuity   

    
    ! Local variables
    integer(ip) :: istat

    ! Parallel environment MUST BE already created
    assert ( p_trian%p_env%created )
    this%p_trian        => p_trian 
	call this%set_dof_descriptor(dof_descriptor)

    if( this%p_trian%p_env%p_context%iam >= 0 ) then
       call serial_fe_space_allocate_structures(  this%serial_fe_space, p_trian%f_trian, dof_descriptor, &
                                                  time_steps_to_store = time_steps_to_store, hierarchical_basis = hierarchical_basis, &
                                                  static_condensation = static_condensation, num_continuity = num_continuity, &
                                                  num_ghosts = p_trian%num_ghosts ) 

       call serial_fe_space_fe_list_create ( this%serial_fe_space, problem, continuity, enable_face_integration, &
            &                                order, material, p_cond%f_conditions )

       ! Communicate problem, continuity, order, and material
       call ghost_elements_exchange ( p_trian%p_env%p_context%icontxt, p_trian%f_el_import, this%serial_fe_space%finite_elements )
       ! Create ghost fem space (only partially, i.e., previous info)
       call ghost_fe_list_create ( this ) 
       call serial_fe_space_integration_faces_list( this%serial_fe_space )
    end if

  end subroutine par_fe_space_create
  
  subroutine par_fe_space_free ( this )
    implicit none
    ! Dummy arguments
    class(par_fe_space_t), intent(inout) :: this  
    
    ! Local variables
    integer(ip) :: ielem
    
    ! Parallel environment MUST BE already created
    assert ( associated(this%p_trian) )
    assert ( this%p_trian%p_env%created )
    
    if( this%p_trian%p_env%p_context%iam >= 0 ) then
       ! Deallocate type(finite_element_ts) associated to ghost elements
       do ielem = this%p_trian%f_trian%num_elems+1, this%p_trian%f_trian%num_elems+this%p_trian%num_ghosts
          if(allocated(this%serial_fe_space%finite_elements(ielem)%reference_element_vars)) call memfree(this%serial_fe_space%finite_elements(ielem)%reference_element_vars,__FILE__,__LINE__)
          if(allocated(this%serial_fe_space%finite_elements(ielem)%elem2dof)) call memfree(this%serial_fe_space%finite_elements(ielem)%elem2dof,__FILE__,__LINE__)
          if(allocated(this%serial_fe_space%finite_elements(ielem)%unkno)) call memfree(this%serial_fe_space%finite_elements(ielem)%unkno,__FILE__,__LINE__)
          call finite_element_free_unpacked(this%serial_fe_space%finite_elements(ielem))
       end do
       call this%serial_fe_space%free()
    end if
    call this%blocks_dof_distribution%free()
    nullify ( this%p_trian )
	nullify ( this%dof_descriptor )
  end subroutine par_fe_space_free

  subroutine par_fe_space_print ( p_fe_space, lunou )
    implicit none
    class(par_fe_space_t), intent(in) :: p_fe_space
	integer(ip)          , intent(in) :: lunou
    if( p_fe_space%p_trian%p_env%p_context%iam >= 0 ) then
	  call p_fe_space%serial_fe_space%print(lunou, p_fe_space%p_trian%num_ghosts)
    end if
  end subroutine par_fe_space_print

  !==================================================================================================
  subroutine par_fe_space_set_analytical_code(p_fe_space,spatial_code,temporal_code)
    implicit none
    class(par_fe_space_t), intent(inout) :: p_fe_space
    integer(ip)         , intent(in)    :: spatial_code(:)
    integer(ip)         , intent(in)    :: temporal_code(:)

    if(p_fe_space%p_trian%p_env%am_i_fine_task()) then
       call p_fe_space%serial_fe_space%set_analytical_code(spatial_code,temporal_code)
    end if

  end subroutine par_fe_space_set_analytical_code

  
  !==================================================================================================
  function par_fe_space_get_blocks_dof_distribution(p_fe_space)
    implicit none
    class(par_fe_space_t), target, intent(in) :: p_fe_space
    type(blocks_dof_distribution_t), pointer :: par_fe_space_get_blocks_dof_distribution 
	par_fe_space_get_blocks_dof_distribution => p_fe_space%blocks_dof_distribution
  end function par_fe_space_get_blocks_dof_distribution

  !*********************************************************************************
  ! This subroutine takes the local finite element space, extended with ghost 
  ! elements that include (after the element exchange) values for order, continuity,
  ! material, problem. With this info, we create the fixed info pointers in ghost
  ! elements, which are needed in order to sort dofs on different processors the 
  ! same way.
  !*********************************************************************************
  subroutine ghost_fe_list_create( p_fe_space )
    implicit none
    type(par_fe_space_t), target, intent(inout) :: p_fe_space

    integer(ip) :: ielem, nvars, f_type, ivar, f_order, istat, pos_elinf, v_key
    logical :: created
    integer(ip) :: aux_val, max_num_nodes, lndof, nnode

    do ielem = p_fe_space%p_trian%f_trian%num_elems+1, p_fe_space%p_trian%f_trian%num_elems+p_fe_space%p_trian%num_ghosts
       !write (*,*) '************* GHOST ELEMENT *************',ielem
       nvars = p_fe_space%serial_fe_space%dof_descriptor%problems(p_fe_space%serial_fe_space%finite_elements(ielem)%problem)%p%nvars
       p_fe_space%serial_fe_space%finite_elements(ielem)%num_vars = nvars
       f_type = p_fe_space%p_trian%f_trian%elems(ielem)%geo_reference_element%ftype
       !write(*,*) 'f_type ghosts',f_type
       !write(*,*) 'nvars',nvars
       assert ( f_type > 0)
       call memalloc(nvars, p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars, __FILE__, __LINE__ )
       allocate(p_fe_space%serial_fe_space%finite_elements(ielem)%nodes_per_vef(nvars), stat=istat )
       do ivar=1,nvars
          f_order = p_fe_space%serial_fe_space%finite_elements(ielem)%order(ivar)
          !write(*,*) 'f_order',f_order
          v_key = p_fe_space%p_trian%f_trian%num_dims + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
          call p_fe_space%serial_fe_space%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          if ( istat == new_index) then 
             !write (*,*) ' FIXED INFO NEW'
             call reference_element_create(p_fe_space%serial_fe_space%finite_elements_info(pos_elinf),f_type,              &
                  &                             f_order,p_fe_space%p_trian%f_trian%num_dims)
          end if
          p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p => p_fe_space%serial_fe_space%finite_elements_info(pos_elinf)
          if ( p_fe_space%serial_fe_space%finite_elements(ielem)%continuity(ivar) /= 0 ) then
             p_fe_space%serial_fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p => p_fe_space%serial_fe_space%finite_elements_info(pos_elinf)%ndxob
          else 
             p_fe_space%serial_fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p => p_fe_space%serial_fe_space%finite_elements_info(pos_elinf)%ndxob_int
          end if
       end do

       ! Compute number of DOFs in the elemental matrix associated to ielem
       lndof = 0
       max_num_nodes = 0
       do ivar=1,nvars
          if ( p_fe_space%serial_fe_space%static_condensation ) then
             nnode = p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode -                                   &
                  &  p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nodes_vef(p_fe_space%serial_fe_space%g_trian%num_dims+1) ! SB.alert : do not use nodes_vef
          else
             nnode = p_fe_space%serial_fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode 
          end if
          lndof = lndof + nnode
          max_num_nodes = max(max_num_nodes,nnode)
       end do

       call memalloc( max_num_nodes, nvars, p_fe_space%serial_fe_space%finite_elements(ielem)%elem2dof, __FILE__,__LINE__ )
       p_fe_space%serial_fe_space%finite_elements(ielem)%elem2dof = 0
       call memalloc( max_num_nodes, nvars, p_fe_space%serial_fe_space%time_steps_to_store, p_fe_space%serial_fe_space%finite_elements(ielem)%unkno, __FILE__,__LINE__)
       p_fe_space%serial_fe_space%finite_elements(ielem)%unkno = 0.0_rp
    end do

  end subroutine ghost_fe_list_create

  subroutine interface_faces_list( p_trian, p_fe_space ) 
    implicit none
    ! Parameters
    type(par_triangulation_t), intent(in)     :: p_trian 
    type(par_fe_space_t)    , intent(inout)  :: p_fe_space

    integer(ip) :: count_int, iobje, ielem, istat

    ! interface faces (among subdomains)

    count_int = 0
    do iobje = 1, p_trian%f_trian%num_vefs
       if ( p_trian%f_trian%vefs(iobje)%dimension == p_trian%f_trian%num_dims-1 ) then
          if (p_trian%vefs(iobje)%interface == 1 ) then
             assert( p_trian%f_trian%vefs(iobje)%num_elems_around == 2 )
             ielem = p_trian%f_trian%vefs(iobje)%elems_around(1)
             count_int = count_int + 1
          end if
       end if
    end do

    allocate( p_fe_space%interface_faces(count_int), stat=istat)
    check ( istat == 0 )

    count_int = 0
    do iobje = 1, p_trian%f_trian%num_vefs
       if ( p_trian%f_trian%vefs(iobje)%dimension == p_trian%f_trian%num_dims-1 ) then
          if (p_trian%vefs(iobje)%interface == 1 ) then
             assert( p_trian%f_trian%vefs(iobje)%num_elems_around == 2 )
             ielem = p_trian%f_trian%vefs(iobje)%elems_around(1)
             count_int = count_int + 1
             p_fe_space%interface_faces(count_int)%face_vef = iobje
          end if
       end if
    end do

   p_fe_space%num_interface_faces = count_int

  end subroutine interface_faces_list

end module par_fe_space_names

