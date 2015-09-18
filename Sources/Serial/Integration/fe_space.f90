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
  ! Modules
  use types_names
  use memor_names
  use array_names
  use triangulation_names
  use hash_table_names
  use problem_names
  use integration_tools_names
  use interpolation_tools_names
  !use face_integration_names
  use fe_space_types_names
  use dof_descriptor_names
  use migratory_element_names
  use conditions_names
  use finite_element_names
  use graph_names
  use serial_scalar_matrix_names
  use serial_block_matrix_names
  use sort_names


#ifdef memcheck
  use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private

  integer(ip), parameter       :: max_global_interpolations  = 50, &  ! Maximum number of interpolations
                                & max_number_problems  = 2               ! Maximum number of problems

  ! Global information of the fe_space
  type fe_space_t  
     integer(ip)                       :: num_continuity       ! Number of materials (maximum value)
     logical                           :: static_condensation  ! Flag for static condensation 
     logical                           :: hierarchical_basis   ! Flag for hierarchical basis
     class(migratory_element_t), allocatable :: mig_elems(:)         ! Migratory elements list_t
     type(finite_element_t)    , pointer     :: finite_elements(:)             ! List of FEs
     type(fe_face_t)           , allocatable :: fe_faces(:)             ! List of active faces

     type(triangulation_t)  , pointer :: g_trian => NULL() ! Triangulation
     type(dof_descriptor_t)        , pointer :: dof_descriptor

     ! Array of working arrays (element matrix/vector) (to be pointed from finite_elements)
     type(position_hash_table_t)          :: pos_elmatvec
     type(array_rp2_t)                    :: lelmat(max_global_interpolations)
     type(array_rp1_t)                    :: lelvec(max_global_interpolations)

     ! Integrator
     type(position_hash_table_t)          :: pos_volume_integrator
     type(position_hash_table_t)          :: pos_face_integrator
     type(volume_integrator_t)            :: lvoli(max_global_interpolations)        
     type(face_integrator_t)              :: lfaci(max_global_interpolations)

     ! Interpolator
     type(array_rp2_t)                    :: linter(max_global_interpolations)

     ! Starting DOF position
     type(position_hash_table_t)          :: pos_start
     type(array_ip1_t)                    :: lstart(max_number_problems)

     ! Array of reference elements (to be pointed from finite_elements)
     type(position_hash_table_t)          :: pos_elem_info
     type(reference_element_t)            :: finite_elements_info(max_elinf)

     ! Acceleration arrays
     type(list_2d_t), allocatable       :: vef2dof(:)       ! An auxiliary array to accelerate some parts of the code
     integer(ip), allocatable           :: ndofs(:)            ! number of dofs (nblocks)
     integer(ip)                        :: time_steps_to_store ! Time steps to store in unkno

     ! List of faces where we want to integrate things
     type(fe_face_t), allocatable        :: interior_faces(:), boundary_faces(:)
     integer(ip)                         :: num_interior_faces, num_boundary_faces

     ! Analytical function auxiliar array
     type(array_ip2_t)                    :: l_analytical_code(max_number_problems)

     ! Much better here rather than as a module variable
     !type(list_t) :: void_list_t

   contains
     procedure, private :: fe_space_create_make_serial_scalar_coefficient_matrix
     procedure, private :: fe_space_create_make_serial_block_coefficient_matrix
	 generic :: make_coefficient_matrix => fe_space_create_make_serial_scalar_coefficient_matrix, &
	                                       fe_space_create_make_serial_block_coefficient_matrix
     procedure :: set_analytical_code => fe_space_set_analytical_code
  end type fe_space_t

  ! Types
  public :: fe_space_t

  ! Methods
  public :: fe_space_create, fe_space_print, fe_space_free, &
       &    fe_space_integration_faces_list, fe_space_allocate_structures, &
       &    fe_space_fe_list_create, setup_dof_graph_from_block_row_col_identifiers

contains

  subroutine fe_space_create_make_serial_scalar_coefficient_matrix(fe_space,symmetric_storage,is_symmetric,sign,serial_scalar_matrix)
    implicit none
	class(fe_space_t)            , intent(in)  :: fe_space
    logical                      , intent(in)  :: symmetric_storage
    logical                      , intent(in)  :: is_symmetric
	integer(ip)                  , intent(in)  :: sign
    type(serial_scalar_matrix_t) , intent(out) :: serial_scalar_matrix
	
	assert ( fe_space%dof_descriptor%nblocks == 1 ) 
	call serial_scalar_matrix%create(symmetric_storage,is_symmetric,sign)
	call setup_dof_graph_from_block_row_col_identifiers ( 1, 1, fe_space, serial_scalar_matrix%graph )
	call serial_scalar_matrix%allocate()
  end subroutine fe_space_create_make_serial_scalar_coefficient_matrix
  
  subroutine fe_space_create_make_serial_block_coefficient_matrix(fe_space, & 
												  diagonal_blocks_symmetric_storage, & 
												  diagonal_blocks_symmetric, &
												  diagonal_blocks_sign,& 
												  serial_block_matrix)
    implicit none
    class(fe_space_t)            , intent(in)  :: fe_space
	logical                     , intent(in)  :: diagonal_blocks_symmetric_storage(fe_space%dof_descriptor%nblocks)
    logical                     , intent(in)  :: diagonal_blocks_symmetric(fe_space%dof_descriptor%nblocks)
    integer(ip)                 , intent(in)  :: diagonal_blocks_sign(fe_space%dof_descriptor%nblocks)
	type(serial_block_matrix_t) , intent(out) :: serial_block_matrix
	
	
  end subroutine fe_space_create_make_serial_block_coefficient_matrix


  !==================================================================================================
  ! Allocation of variables in fe_space according to the values in g_trian
  subroutine fe_space_create( g_trian, dof_descriptor, fe_space, problem, bcond, continuity, order, & 
                              material, time_steps_to_store, hierarchical_basis,  & 
                              static_condensation, num_continuity, num_ghosts )
    implicit none
    type(triangulation_t), target, intent(in)    :: g_trian   
    type(dof_descriptor_t)      , target, intent(in)    :: dof_descriptor
    type(fe_space_t)        , target, intent(inout) :: fe_space
    integer(ip)                    , intent(in)    :: problem(:)
    type(conditions_t)           , intent(in)    :: bcond
    integer(ip)                    , intent(in)    :: continuity(:,:)
    integer(ip)                    , intent(in)    :: order(:,:)
    integer(ip)                    , intent(in)    :: material(:)
    integer(ip)          , optional, intent(in)    :: time_steps_to_store
    logical          , optional, intent(in)    :: hierarchical_basis
    logical          , optional, intent(in)    :: static_condensation
    integer(ip)          , optional, intent(in)    :: num_continuity   
    integer(ip)          , optional, intent(in)    :: num_ghosts         

    integer(ip) :: istat, num_ghosts_

    call fe_space_allocate_structures( g_trian, dof_descriptor, fe_space, time_steps_to_store = time_steps_to_store,&
         hierarchical_basis = hierarchical_basis, static_condensation = static_condensation, &
          num_continuity = num_continuity, num_ghosts = num_ghosts )  

    call fe_space_fe_list_create ( fe_space, problem, continuity, order, material, bcond )

    call fe_space_integration_faces_list( fe_space )

  end subroutine fe_space_create

  !==================================================================================================
  ! Allocation of variables in fe_space according to the values in g_trian
  subroutine fe_space_allocate_structures( g_trian, dof_descriptor, fe_space, time_steps_to_store, &
       & hierarchical_basis, static_condensation, num_continuity, num_ghosts )
    implicit none
    type(fe_space_t)   ,  target, intent(inout)                :: fe_space
    type(triangulation_t)   , intent(in), target   :: g_trian   
    type(dof_descriptor_t) , intent(in), target           :: dof_descriptor  
    integer(ip), optional, intent(in) :: time_steps_to_store
    logical, optional, intent(in) :: hierarchical_basis
    logical, optional, intent(in) :: static_condensation
    integer(ip), optional, intent(in) :: num_continuity 
    integer(ip), optional, intent(in) :: num_ghosts            

    integer(ip) :: istat, num_ghosts_

    ! Hierarchical flag
    if (present(hierarchical_basis)) then
       fe_space%hierarchical_basis = hierarchical_basis
    else
       fe_space%hierarchical_basis = .false.
    end if

    ! Static condensation flag
    if (present(static_condensation)) then
       fe_space%static_condensation = static_condensation
    else
       fe_space%static_condensation = .false.
    end if

    ! Number materials flag
    if (present(num_continuity)) then
       fe_space%num_continuity = num_continuity
    else
       fe_space%num_continuity = 1
    end if

    ! Time steps to store flag
    if (present(time_steps_to_store)) then
       fe_space%time_steps_to_store = time_steps_to_store
    else
       fe_space%time_steps_to_store = 1
    end if

    ! Ghosts elements (parallel case)
    if (present(num_ghosts)) then
       num_ghosts_ = num_ghosts
    else
       num_ghosts_ = 0
    end if

    allocate( finite_element_t :: fe_space%mig_elems(g_trian%num_elems + num_ghosts_), stat=istat)
    check ( istat == 0 )
    
    select type( this => fe_space%mig_elems )
    type is(finite_element_t)
       fe_space%finite_elements => this
    end select

    !  Initialization of pointer to triangulation
    fe_space%g_trian => g_trian

    !  Initialization of pointer to dof_descriptor
    fe_space%dof_descriptor => dof_descriptor

    ! Allocation of elemental matrix and vector parameters
    call fe_space%pos_elmatvec%init(ht_length)

    ! Initialization of volume and face integrators parameters
    call fe_space%pos_volume_integrator%init(ht_length)
    call fe_space%pos_face_integrator%init(ht_length)

    ! Initialization of element fixed info parameters
    call fe_space%pos_elem_info%init(ht_length)

    ! Initialization of starting DOF position array
    call fe_space%pos_start%init(ht_length)

  end subroutine fe_space_allocate_structures

  !==================================================================================================
  ! Fill the fe_space_t assuming that all elements are of type f_type but each variable has different
  ! interpolation order
  subroutine fe_space_fe_list_create( fe_space, problem, continuity, order, material, bcond )
    implicit none
    type(fe_space_t), intent(inout), target  :: fe_space
    integer(ip)    , intent(in)       :: material(:), order(:,:), problem(:)
    integer(ip)    , intent(in)       :: continuity(:,:)
    type(conditions_t), intent(in)  :: bcond

    integer(ip) :: nunk, v_key, ltype(2), nnode, max_num_nodes, nunk_tot, dim, f_order, f_type, nvars, nvars_tot
    integer(ip) :: ielem, istat, pos_elmatvec, pos_elinf, pos_voint, ivar, lndof, iobje
    integer(ip) :: aux_val

    ! Loop over elements
    nvars_tot = fe_space%dof_descriptor%nvars_global
    dim = fe_space%g_trian%num_dims

    ! Material
    fe_space%finite_elements(1:fe_space%g_trian%num_elems)%material = material

    ! Continuity
    do ielem = 1, fe_space%g_trian%num_elems
       ! Assign type of problem and approximation to ielem
       fe_space%finite_elements(ielem)%problem =  problem(ielem)

       nvars = fe_space%dof_descriptor%problems(problem(ielem))%p%nvars
       fe_space%finite_elements(ielem)%num_vars = nvars
       f_type = fe_space%g_trian%elems(ielem)%geo_reference_element%ftype

       ! Set continuity per unknown
       call memalloc(nvars, fe_space%finite_elements(ielem)%continuity, __FILE__, __LINE__)

       ! Set order per unknown
       call memalloc(nvars, fe_space%finite_elements(ielem)%order, __FILE__, __LINE__)

       ! Allocate the fixed info array
       call memalloc(nvars, fe_space%finite_elements(ielem)%reference_element_vars, __FILE__, __LINE__)

       ! Allocate nodes per vef info
       allocate(fe_space%finite_elements(ielem)%nodes_per_vef(nvars), stat=istat )
       check( istat==0 )

       ! Assign pointer to interpolation fixed information and nodes per vef
       do ivar=1,nvars

          ! JP: indices of these arrays (continuity and order) should be changed to (nvars,nelem)
          fe_space%finite_elements(ielem)%continuity(ivar) = continuity(ielem,fe_space%dof_descriptor%problems(problem(ielem))%p%l2g_var(ivar))
          
          fe_space%finite_elements(ielem)%order(ivar) = order(ielem,fe_space%dof_descriptor%problems(problem(ielem))%p%l2g_var(ivar))
          if ( fe_space%finite_elements(ielem)%continuity(ivar) /= 0 ) then
             fe_space%static_condensation = .false. ! Static condensation + dG not possible
          end if

          f_order = fe_space%finite_elements(ielem)%order(ivar)
          v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
          call fe_space%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          if ( istat == new_index) then 
             call reference_element_create(fe_space%finite_elements_info(pos_elinf),f_type,              &
                  &                             f_order,dim)
          end if
          fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p => fe_space%finite_elements_info(pos_elinf)

          if ( fe_space%finite_elements(ielem)%continuity(ivar) /= 0 ) then
             fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p => fe_space%finite_elements_info(pos_elinf)%ndxob
          else 
             fe_space%finite_elements(ielem)%nodes_per_vef(ivar)%p => fe_space%finite_elements_info(pos_elinf)%ndxob_int
          end if
       end do

       ! Assign pointer to geometrical fixed information (assumed to be of order 1)
       v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       call fe_space%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
       if ( istat == new_index) then 
          call reference_element_create(fe_space%finite_elements_info(pos_elinf),f_type,              &
               &                             1,dim)
       end if
       fe_space%finite_elements(ielem)%p_geo_reference_element => fe_space%finite_elements_info(pos_elinf)

       ! Compute number of DOFs in the elemental matrix associated to ielem
       lndof = 0
       max_num_nodes = 0
       do ivar=1,nvars
          if ( fe_space%static_condensation ) then
             ! JP: Is this working? Because max_num_nodes is sent to integ_create to build the interpolation
             !     which is needed even with static_condensation...
             ! JP & SB: In fact we also need to think about it for the unknown. Only elem2dof can be reduced
             !          to the interface.
             nnode = fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode -                                   &
                  &  fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nodes_vef(dim+1) ! SB.alert : do not use nodes_vef
          else
             nnode = fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode 
          end if
          lndof = lndof + nnode
          max_num_nodes = max(max_num_nodes,nnode)
       end do

       ! Assign pointer to p_mat and p_vec in ielem
       call fe_space%pos_elmatvec%get(key=lndof,val=pos_elmatvec,stat=istat)
       if ( istat == new_index ) then
          call array_create ( lndof, lndof, fe_space%lelmat(pos_elmatvec) )
          call array_create ( lndof, fe_space%lelvec(pos_elmatvec) )
       end if
       fe_space%finite_elements(ielem)%p_mat => fe_space%lelmat(pos_elmatvec)
       fe_space%finite_elements(ielem)%p_vec => fe_space%lelvec(pos_elmatvec)

       ! Allocate elem2dof, unkno, bc_code
       call memalloc( max_num_nodes, nvars, fe_space%finite_elements(ielem)%elem2dof, __FILE__,__LINE__ )
       fe_space%finite_elements(ielem)%elem2dof = 0
       call memalloc( max_num_nodes, nvars, fe_space%time_steps_to_store, fe_space%finite_elements(ielem)%unkno, __FILE__,__LINE__)
       fe_space%finite_elements(ielem)%unkno = 0.0_rp
       call memalloc(nvars,fe_space%finite_elements(ielem)%integ,__FILE__,__LINE__)
       call memalloc(nvars,fe_space%finite_elements(ielem)%inter,__FILE__,__LINE__)
       call memalloc(nvars,fe_space%finite_elements(ielem)%p_geo_reference_element%nvef,fe_space%finite_elements(ielem)%bc_code,__FILE__,__LINE__, 0)

       ! Fill bc_code
       do iobje = 1,fe_space%finite_elements(ielem)%p_geo_reference_element%nvef
          do ivar=1,nvars
             fe_space%finite_elements(ielem)%bc_code(ivar,iobje) = bcond%code( fe_space%dof_descriptor%problems(problem(ielem))%p%l2g_var(ivar),  fe_space%g_trian%elems(ielem)%vefs(iobje) )
          end do
       end do

       ! Assign pointers to volume integration & interpolator
       ltype(2) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       do ivar = 1,nvars
          f_order  = fe_space%finite_elements(ielem)%order(ivar)
          ltype(1) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
          v_key    = (max_ndime+1)*(max_FE_types+1)*(max_order) * ltype(1) + ltype(2)
          call fe_space%pos_volume_integrator%get(key=v_key, val=pos_voint, stat = istat)
          ! SB.alert : g_ord = 1 !!!! But only for linear geometry representation
          if ( istat == new_index ) then
             call volume_integrator_create(f_type,f_type,dim,1,f_order,fe_space%lvoli(pos_voint),     &
                  &                        khie = fe_space%hierarchical_basis, mnode=max_num_nodes)
             call interpolator_create(f_type,f_type,dim,1,f_order,                                          &
                  &                   fe_space%finite_elements(ielem)%p_geo_reference_element%nnode,        &
                  &                   fe_space%finite_elements(ielem)%reference_element_vars(ivar)%p%nnode, &
                  &                   fe_space%linter(pos_voint),khie = fe_space%hierarchical_basis)
          end if
          fe_space%finite_elements(ielem)%integ(ivar)%p => fe_space%lvoli(pos_voint)
          fe_space%finite_elements(ielem)%inter(ivar)%p => fe_space%linter(pos_voint)
       end do

       ! Assign pointers to starting DOF position & analytical code
       call fe_space%pos_start%get(key=nvars, val=pos_voint, stat = istat)
       if ( istat == new_index ) then
          call pointer_variable(fe_space%finite_elements(ielem),fe_space%dof_descriptor,fe_space%lstart(pos_voint))
          call array_create(nvars,2,fe_space%l_analytical_code(pos_voint))
       end if
       fe_space%finite_elements(ielem)%start             => fe_space%lstart(pos_voint)
       fe_space%finite_elements(ielem)%p_analytical_code => fe_space%l_analytical_code(pos_voint)
       fe_space%l_analytical_code(pos_voint)%a = 0
    end do

  end subroutine fe_space_fe_list_create

  !==================================================================================================
  subroutine dof_to_elnod( el2dof, idof, ivar, nvef, pnodob, lnode )
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: el2dof(:,:)
    integer(ip), intent(in)  :: idof, ivar, nvef
    integer(ip), intent(in)  :: pnodob(nvef+1)
    integer(ip), intent(out) :: lnode

    ! Locals
    integer(ip) :: iobje, jnode

    do iobje = 1,nvef
       do jnode = pnodob(iobje),pnodob(iobje+1)-1
          if ( el2dof(jnode,ivar) == idof ) then
             lnode = jnode - pnodob(iobje) + 1
             return
          end if
       end do
    end do

    ! ERROR: not found !!!!
    assert(0==1)

  end subroutine dof_to_elnod

  subroutine fe_space_print ( lunou, fe_space, num_ghosts )
    implicit none
    integer(ip)      , intent(in) :: lunou
    type(fe_space_t), intent(in) :: fe_space
    integer(ip), optional, intent(in) :: num_ghosts    

    integer(ip) :: ielem, num_ghosts_
    
    ! Ghosts elements (parallel case)
    if (present(num_ghosts)) then
       num_ghosts_ = num_ghosts
    else
       num_ghosts_ = 0
    end if
    
    write (lunou,*) 'Number of materials: ', fe_space%num_continuity
    write (lunou,*) 'Static condensation flag: ', fe_space%static_condensation

    write (lunou,*) '****PRINT ELEMENT LIST INFO****'
    do ielem = 1, fe_space%g_trian%num_elems + num_ghosts_
       write (lunou,*) '****PRINT ELEMENT ',ielem,' INFO****'
       call finite_element_print ( lunou, fe_space%finite_elements(ielem) )
       write (lunou,*) '****END PRINT ELEMENT ',ielem,' INFO****'
    end do

    write (lunou,*) 'Number boundary faces: ', fe_space%num_boundary_faces
    write (lunou,*) 'Number interior faces: ', fe_space%num_interior_faces
    write (lunou,*) 'Boundary faces: ', fe_space%boundary_faces(:)%face_vef
    write (lunou,*) 'Interior faces: ', fe_space%interior_faces(:)%face_vef
    

  end subroutine fe_space_print

 !==================================================================================================
  subroutine fe_space_free ( fe_space )
    implicit none
    type(fe_space_t), intent(inout) :: fe_space
    integer(ip)                    :: i,j
    integer(ip) :: istat

    ! SB.alert : To be thought for the new structure triangulation
    ! lface_free
    ! if (allocated(fe_space%fe_faces)) then
    !    j = 1
    !    do i = 1, fe_space%g_trian%num_vefs
    !       if ( trian%vefs(i)%dimension == fe_space%g_trian%num_dims-1 ) then
    !          fe_space%fe_faces(i)%nei_elem  = 0
    !          fe_space%fe_faces(i)%pos_elem  = 0
    !          fe_space%fe_faces(i)%subface   = 0
    !          fe_space%fe_faces(i)%refinement_level = 0
    !          nullify (fe_space%fe_faces(i)%p_mat)
    !          nullify (fe_space%fe_faces(i)%p_vec)
    !          if (allocated(fe_space%fe_faces(i)%integ)) call memfree(fe_space%fe_faces(i)%integ,__FILE__,__LINE__)
    !          if ( trian%vefs(i)%border /= -1) then
    !             do j=1,fe_space%finite_elements(fe_space%fe_faces(i)%nei_elem(1))%p%nvars
    !                call array_free(fe_space%fe_faces(i)%o2n(j))
    !             end do
    !             call memfree(fe_space%fe_faces(i)%o2n,__FILE__,__LINE__)
    !          end if
    !          j = j+1
    !       end if
    !    end do

    !    deallocate(fe_space%fe_faces)
    !    ! Face integrators free
    !    do i = 1,fe_space%cur_lfaci-1 
    !       call integ_free( fe_space%lfaci(i) )
    !    end do
    !    fe_space%cur_lfaci = 0
    !    call fe_space%ht_pos_face_integrator%free
    ! end if

    

    do i = 1, fe_space%g_trian%num_elems
!       nullify ( fe_space%finite_elements(i)%problem )
       call memfree( fe_space%finite_elements(i)%elem2dof,__FILE__,__LINE__)
       call memfree( fe_space%finite_elements(i)%unkno,__FILE__,__LINE__)
       !call memfree( fe_space%finite_elements(i)%jvars,__FILE__,__LINE__)
       if (allocated(fe_space%finite_elements(i)%bc_code)) then
          call memfree(fe_space%finite_elements(i)%bc_code,__FILE__,__LINE__)
       end if
       nullify ( fe_space%finite_elements(i)%p_geo_reference_element )
       nullify ( fe_space%finite_elements(i)%p_mat )
       nullify ( fe_space%finite_elements(i)%p_vec )
       nullify ( fe_space%finite_elements(i)%start )
       do j = 1, fe_space%finite_elements(i)%num_vars
          nullify ( fe_space%finite_elements(i)%nodes_per_vef(j)%p )  
       end do
       !if(allocated(fe_space%finite_elements(i)%nodes_per_vef)) call memfree(fe_space%finite_elements(i)%nodes_per_vef,__FILE__,__LINE__)
       if(allocated(fe_space%finite_elements(i)%reference_element_vars)) call memfree(fe_space%finite_elements(i)%reference_element_vars,__FILE__,__LINE__)
       if(allocated(fe_space%finite_elements(i)%integ)) call memfree(fe_space%finite_elements(i)%integ,__FILE__,__LINE__)
       if(allocated(fe_space%finite_elements(i)%inter)) call memfree(fe_space%finite_elements(i)%inter,__FILE__,__LINE__)
       !if(allocated(fe_space%finite_elements(i)%iv))    call memfree(fe_space%finite_elements(i)%iv   ,__FILE__,__LINE__)
       if(allocated(fe_space%finite_elements(i)%continuity))    call memfree(fe_space%finite_elements(i)%continuity   ,__FILE__,__LINE__)
       if(allocated(fe_space%finite_elements(i)%order))    call memfree(fe_space%finite_elements(i)%order   ,__FILE__,__LINE__)
       !if(allocated(fe_space%finite_elements(i)%material))    call memfree(fe_space%finite_elements(i)%iv   ,__FILE__,__LINE__)
       !if(allocated(fe_space%finite_elements(i)%p_nod)) call memfree(fe_space%finite_elements(i)%p_nod,__FILE__,__LINE__)
    end do
    
    deallocate ( fe_space%mig_elems, stat=istat )
    check(istat==0)
    nullify ( fe_space%finite_elements )

    do i = 1,fe_space%pos_elmatvec%last()
       call array_free( fe_space%lelmat(i) )
       call array_free( fe_space%lelvec(i) )
    end do
    !deallocate ( fe_space%lelmat )
    !deallocate ( fe_space%lelvec )
    call fe_space%pos_elmatvec%free
    call fe_space%pos_elmatvec%free

    do i = 1,fe_space%pos_volume_integrator%last()
       call volume_integrator_free( fe_space%lvoli(i) )
       call interpolator_free( fe_space%linter(i) )
    end do
    !deallocate ( fe_space%lvoli )
    call fe_space%pos_volume_integrator%free

    do i = 1,fe_space%pos_elem_info%last()
       call reference_element_free (fe_space%finite_elements_info(i))
    end do
    call fe_space%pos_elem_info%free

    do i = 1,fe_space%pos_start%last()
       call array_free( fe_space%lstart(i) )
       call array_free( fe_space%l_analytical_code(i) )
    end do
    call fe_space%pos_start%free

    nullify ( fe_space%g_trian )

    !call memfree( fe_space%void_list%p, __FILE__, __LINE__ )
    !call memfree( fe_space%void_list%l, __FILE__, __LINE__ )
    !fe_space%void_list%n = 0

    do i = 1, fe_space%dof_descriptor%nblocks
       call memfree (fe_space%vef2dof(i)%p , __FILE__, __LINE__ )
       call memfree (fe_space%vef2dof(i)%l , __FILE__, __LINE__ )
    end do

    call memfree (fe_space%ndofs , __FILE__, __LINE__ )

    deallocate( fe_space%interior_faces, stat=istat)
    check ( istat == 0 )

    deallocate( fe_space%boundary_faces, stat=istat)
    check ( istat == 0 )

  end subroutine fe_space_free

  subroutine get_p_faces ( fe_space, trian )
    implicit none
    type(fe_space_t), intent(in)  :: fe_space
    type(triangulation_t), intent(inout)    :: trian

    !   integer(ip) :: i


    !   num_faces

    !   call memalloc ( trian%num_elems, num_faces, __FILE__,__LINE__ )

    !   do i=1,mesh%nelem
    !         mesh%p_face(i) = mesh%pnods(i) + fe_space%finite_elements(i)%reference_element_vars(1)%p%nvef_dim(mesh%ndime) -1
    !   end do

  end subroutine get_p_faces

  subroutine fe_space_integration_faces_list( fe_space ) 
    implicit none
    ! Parameters
    type(fe_space_t), intent(inout)               :: fe_space

    integer(ip) :: count_int, count_bou, mat_i, mat_j, iobje, ielem, jelem, istat
    integer(ip) :: g_var, iprob, jprob, ivars, jvars

    ! integration faces (interior / boundary)
    ! The list of boundary faces includes all faces, whereas the interior ones are only those
    ! where we expect to integrate things (based on continuity flags)
    count_int = 0
    count_bou = fe_space%g_trian%num_boundary_faces
    do iobje = 1, fe_space%g_trian%num_vefs
       if ( fe_space%g_trian%vefs(iobje)%dimension == fe_space%g_trian%num_dims-1 ) then
          if ( fe_space%g_trian%vefs(iobje)%border == -1 ) then
             assert( fe_space%g_trian%vefs(iobje)%num_elems_around == 2 )
             ielem = fe_space%g_trian%vefs(iobje)%elems_around(1)
             jelem = fe_space%g_trian%vefs(iobje)%elems_around(2)
             iprob = fe_space%finite_elements(ielem)%problem
             jprob = fe_space%finite_elements(jelem)%problem
             do ivars = 1, fe_space%dof_descriptor%problems(iprob)%p%nvars
                g_var = fe_space%dof_descriptor%problems(iprob)%p%l2g_var(ivars)
                jvars = fe_space%dof_descriptor%g2l_vars(g_var,jprob)
                mat_i = fe_space%finite_elements(ielem)%continuity(ivars)
                mat_j = fe_space%finite_elements(jelem)%continuity(jvars)
                if ( mat_i == 0 .or. mat_i /= mat_j ) then
                   count_int = count_int + 1
                   exit
                   !fe_space%interior_faces(count_int) = iobje
                end if
             end do
          else
             assert( fe_space%g_trian%vefs(iobje)%num_elems_around == 1 )
          end if
       end if
    end do

    allocate( fe_space%interior_faces(count_int), stat=istat)
    check ( istat == 0 )
    allocate( fe_space%boundary_faces(count_bou), stat=istat)
    check ( istat == 0 )


    fe_space%num_interior_faces = count_int
    fe_space%num_boundary_faces = count_bou

    count_int = 0
    count_bou = 0
    do iobje = 1, fe_space%g_trian%num_vefs
       if ( fe_space%g_trian%vefs(iobje)%dimension == fe_space%g_trian%num_dims-1 ) then
          if ( fe_space%g_trian%vefs(iobje)%border == -1 ) then
             assert( fe_space%g_trian%vefs(iobje)%num_elems_around == 2 )
             ielem = fe_space%g_trian%vefs(iobje)%elems_around(1)
             jelem = fe_space%g_trian%vefs(iobje)%elems_around(2)
             iprob = fe_space%finite_elements(ielem)%problem
             jprob = fe_space%finite_elements(jelem)%problem
             do ivars = 1, fe_space%dof_descriptor%problems(iprob)%p%nvars
                g_var = fe_space%dof_descriptor%problems(iprob)%p%l2g_var(ivars)
                jvars = fe_space%dof_descriptor%g2l_vars(g_var,jprob)
                mat_i = fe_space%finite_elements(ielem)%continuity(ivars)
                mat_j = fe_space%finite_elements(jelem)%continuity(jvars)
                if ( mat_i == 0 .or. mat_i /= mat_j ) then
                   count_int = count_int + 1
                   fe_space%interior_faces(count_int)%face_vef = iobje
                   exit
                end if
             end do
          else
             assert( fe_space%g_trian%vefs(iobje)%num_elems_around == 1 )
             count_bou = count_bou + 1
             fe_space%boundary_faces(count_bou)%face_vef = iobje
          end if
       end if
    end do


  end subroutine fe_space_integration_faces_list

  !==================================================================================================
  subroutine pointer_variable( finite_element, dof_descriptor, start ) 
    implicit none
	type(finite_element_t), intent(in)  :: finite_element
    type(dof_descriptor_t), intent(in)  :: dof_descriptor
    type(array_ip1_t)     , intent(out) :: start

    integer(ip) :: ivar

    call array_create(dof_descriptor%problems(finite_element%problem)%p%nvars+1,start)

    do ivar = 1,dof_descriptor%problems(finite_element%problem)%p%nvars
       start%a(ivar+1) = finite_element%reference_element_vars(ivar)%p%nnode
    end do
    start%a(1) = 1
    do ivar = 2, dof_descriptor%problems(finite_element%problem)%p%nvars+1
       start%a(ivar) = start%a(ivar) + start%a(ivar-1)
    end do
  end subroutine pointer_variable

  !==================================================================================================
  subroutine fe_space_set_analytical_code(fe_space,spatial_code,temporal_code)
    implicit none
    class(fe_space_t), intent(inout) :: fe_space
    integer(ip)      , intent(in)    :: spatial_code(:)
    integer(ip)      , intent(in)    :: temporal_code(:)
    ! Locals
    integer(ip) :: nvars,position,istat
    
    nvars = size(spatial_code,1)

    ! Checks
    check(nvars==size(temporal_code,1))

    ! Copy codes
    call fe_space%pos_start%get(key=nvars, val=position, stat = istat)
    if ( istat == old_index ) then
       fe_space%l_analytical_code(position)%a(:,1) = spatial_code
       fe_space%l_analytical_code(position)%a(:,2) = temporal_code
    else
       write(*,*) 'fe_space%set_analytical_code: unsupported size of spatial codes array'
       check(.false.)
    end if

  end subroutine fe_space_set_analytical_code
  
  !*********************************************************************************
  ! This subroutine takes the fe_space and creates a dof_graph. 
  ! The dof_graph includes both the coupling by continuity like in continuous Galerkin 
  ! methods, and the coupling by face terms (of discontinuous Galerkin type). The algorithm 
  ! considers both the case with static condensation and without it. In order to call this 
  ! subroutine, we need to compute first element2dof and vef2dof arrays.
  ! 3) It generates the local graph (to be extended in *par_create_global_dof_info_names*
  !    to put additional DOFs due to face integration coupling with DOFs from ghost 
  !    elements) (see explanation of the subroutine below and *ghost_dofs_by_integration*
  !    in *par_create_global_dof_info_names* for more insight)
  !*********************************************************************************
  subroutine setup_dof_graph_from_block_row_col_identifiers( iblock, jblock, fe_space, dof_graph ) 
    implicit none
    ! Parameters
    integer(ip)            , intent(in)     :: iblock, jblock
    type(fe_space_t)       , intent(in)     :: fe_space 
    type(graph_t)          , intent(inout)  :: dof_graph

    ! Local variables
    integer(ip) :: iprob, l_var, count, iobje, ielem, jelem, nvapb, ivars, g_var, inter, inode, l_node
    integer(ip) :: idof, jdof, int_i, int_j, istat, jnode, job_g, jobje, jvars, k_var , touch
    integer(ip) :: l_dof, m_dof, m_node, m_var, posi, posf, l_mat, m_mat, knode

    integer(ip) :: nvapbi, nvapbj, nnode, i, iface, jprob, l_faci, l_facj, ic


    integer(ip), allocatable :: aux_ia(:)
    type(hash_table_ip_ip_t) :: visited
	
    touch = 1

    ! Initialize
    dof_graph%nv  = fe_space%ndofs(iblock)
    dof_graph%nv2 = fe_space%ndofs(jblock)
    call memalloc( dof_graph%nv+1, dof_graph%ia, __FILE__,__LINE__ )
    dof_graph%ia = 0

    ! COUNT PART
    call count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, fe_space%dof_descriptor, fe_space%g_trian, fe_space, dof_graph )
    call count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, fe_space%dof_descriptor, fe_space%g_trian, fe_space, dof_graph ) 
    call count_nnz_all_dofs_vs_all_dofs_by_face_integration( iblock, jblock, fe_space%dof_descriptor, fe_space%g_trian, fe_space, dof_graph )  

    !
    dof_graph%ia(1) = 1
    do idof = 2, fe_space%ndofs(iblock)+1
       dof_graph%ia(idof) = dof_graph%ia(idof) + dof_graph%ia(idof-1)
    end do

    call memalloc ( dof_graph%ia(fe_space%ndofs(iblock)+1)-1, dof_graph%ja, __FILE__, __LINE__ )

    ! LIST PART
    call memalloc( dof_graph%nv+1, aux_ia, __FILE__,__LINE__ )
    aux_ia = dof_graph%ia

    call list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, fe_space%dof_descriptor, fe_space%g_trian, fe_space, dof_graph, aux_ia ) 
    call list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, fe_space%dof_descriptor, fe_space%g_trian, fe_space, dof_graph, aux_ia )
    call list_nnz_all_dofs_vs_all_dofs_by_face_integration( iblock, jblock, fe_space%dof_descriptor, fe_space%g_trian, fe_space, dof_graph, aux_ia ) 


    do idof = 1, fe_space%ndofs(iblock)
       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       posi = dof_graph%ia(idof)
       posf = dof_graph%ia(idof+1)-1
       call sort(posf-posi+1,dof_graph%ja(posi:posf))
    end do

    call memfree (aux_ia,__FILE__,__LINE__)

  end subroutine setup_dof_graph_from_block_row_col_identifiers

  !*********************************************************************************
  ! Count NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
  ! both interior and interface nodes.
  !*********************************************************************************
  subroutine count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dof_descriptor, trian, fe_space, dof_graph )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_descriptor_t), intent(in)               :: dof_descriptor
    type(triangulation_t), intent(in)         :: trian 
    type(fe_space_t), intent(in)                 :: fe_space 
    type(graph_t), intent(inout)                :: dof_graph

    ! Local variables
    type(hash_table_ip_ip_t) :: visited
    integer(ip) :: idof, ielem, inode, iobje, iprob, istat, ivars 
    integer(ip) :: jdof, jelem, job_g, jobje, k_var, l_dof, l_mat
    integer(ip) :: l_node, l_var, m_dof, m_mat, m_var, nvapb, touch

    do iobje = 1, trian%num_vefs             
       if ( fe_space%vef2dof(iblock)%p(iobje+1)-fe_space%vef2dof(iblock)%p(iobje) > 0) then
          call visited%init(100) 
          do ielem = 1, trian%vefs(iobje)%num_elems_around
             jelem = trian%vefs(iobje)%elems_around(ielem)
             if ( jelem <= trian%num_elems ) then
                do jobje = 1, trian%elems(jelem)%num_vefs
                   job_g = trian%elems(jelem)%vefs(jobje)
                   call visited%put(key=job_g, val=touch, stat=istat)
                   if ( istat == now_stored ) then   ! interface-interface
                      do idof = fe_space%vef2dof(iblock)%p(iobje), fe_space%vef2dof(iblock)%p(iobje+1)-1
                         l_dof = fe_space%vef2dof(iblock)%l(idof,1)
                         l_var = fe_space%vef2dof(iblock)%l(idof,2)
                         l_mat = fe_space%vef2dof(iblock)%l(idof,3)
                         do jdof = fe_space%vef2dof(jblock)%p(job_g), fe_space%vef2dof(jblock)%p(job_g+1)-1
                            m_dof = fe_space%vef2dof(jblock)%l(jdof,1)
                            m_var = fe_space%vef2dof(jblock)%l(jdof,2)
                            m_mat = fe_space%vef2dof(jblock)%l(jdof,3)
                            if ( dof_descriptor%dof_coupl(l_var,m_var) == 1 .and. l_mat == m_mat ) then
                               if ( .not. dof_graph%symmetric_storage ) then
                                  dof_graph%ia(l_dof+1) = &
                                       & dof_graph%ia(l_dof+1) + 1
                               else
                                  if ( m_dof >= l_dof ) then
                                     dof_graph%ia(l_dof+1) = &
                                          & dof_graph%ia(l_dof+1) + 1
                                  end if
                               end if
                            end if
                         end do
                      end do
                   end if
                end do
                !end do
                if (.not.fe_space%static_condensation) then  ! interface-interior
                   iprob = fe_space%finite_elements(jelem)%problem
                   nvapb = dof_descriptor%prob_block(jblock,iprob)%nd1
                   do idof = fe_space%vef2dof(iblock)%p(iobje), fe_space%vef2dof(iblock)%p(iobje+1)-1
                      l_dof = fe_space%vef2dof(iblock)%l(idof,1)
                      l_var = fe_space%vef2dof(iblock)%l(idof,2)
                      l_mat = fe_space%vef2dof(iblock)%l(idof,3)
                      do ivars = 1, nvapb
                         k_var = dof_descriptor%prob_block(jblock,iprob)%a(ivars)
                         m_var = dof_descriptor%problems(iprob)%p%l2g_var(k_var)
                         m_mat = fe_space%finite_elements(jelem)%continuity(m_var)
                         if ( dof_descriptor%dof_coupl(l_var, m_var) == 1 .and. l_mat == m_mat ) then                
                            if ( .not. dof_graph%symmetric_storage ) then
                               dof_graph%ia(l_dof+1) =  dof_graph%ia(l_dof+1) &
                                    & + fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(jobje+1) &
                                    & - fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(jobje)
                            else 
                               do inode = fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(jobje), &
                                    & fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(jobje+1)-1
                                  l_node = fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%l(inode)
                                  m_dof = fe_space%finite_elements(jelem)%elem2dof(l_node,k_var)
                                  if ( m_dof >= l_dof ) then
                                     dof_graph%ia(l_dof+1) = &
                                          & dof_graph%ia(l_dof+1) + 1
                                  end if
                               end do
                            end if
                         end if
                      end do
                   end do
                end if
             end if
          end do
          call visited%free
       end if
    end do
  end subroutine count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity

  !*********************************************************************************
  ! List NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
  ! both interior and interface nodes.
  !*********************************************************************************
  subroutine list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dof_descriptor, trian, fe_space, dof_graph, aux_ia )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_descriptor_t), intent(in)               :: dof_descriptor
    type(triangulation_t), intent(in)         :: trian 
    type(fe_space_t), intent(in)                 :: fe_space 
    type(graph_t), intent(inout)                :: dof_graph
    integer(ip), intent(inout)                :: aux_ia(:)

    ! Local variables
    type(hash_table_ip_ip_t) :: visited
    integer(ip) :: idof, ielem, inode, iobje, iprob, istat, ivars 
    integer(ip) :: jdof, jelem, job_g, jobje, k_var, l_dof, l_mat
    integer(ip) :: l_node, l_var, m_dof, m_mat, m_var, nvapb, touch, count, ic

    count = 0
    do iobje = 1, trian%num_vefs 
       if ( fe_space%vef2dof(iblock)%p(iobje+1)-fe_space%vef2dof(iblock)%p(iobje) > 0) then
          call visited%init(100) 
          do ielem = 1, trian%vefs(iobje)%num_elems_around
             jelem = trian%vefs(iobje)%elems_around(ielem)
             if ( jelem <= trian%num_elems ) then 
                do jobje = 1, trian%elems(jelem)%num_vefs
                   job_g = trian%elems(jelem)%vefs(jobje)
                   call visited%put(key=job_g, val=touch, stat=istat)
                   if ( istat == now_stored ) then  ! interface-interface
                      do idof = fe_space%vef2dof(iblock)%p(iobje), fe_space%vef2dof(iblock)%p(iobje+1)-1
                         l_dof = fe_space%vef2dof(iblock)%l(idof,1)
                         l_var = fe_space%vef2dof(iblock)%l(idof,2)
                         do jdof = fe_space%vef2dof(jblock)%p(job_g), fe_space%vef2dof(jblock)%p(job_g+1)-1
                            m_dof = fe_space%vef2dof(jblock)%l(jdof,1)
                            m_var = fe_space%vef2dof(jblock)%l(jdof,2)
                            if ( dof_descriptor%dof_coupl(l_var,m_var) == 1 ) then
                               if ( .not. dof_graph%symmetric_storage ) then
                                  ic = aux_ia(l_dof)
                                  dof_graph%ja(ic) = m_dof
                                  aux_ia(l_dof) = aux_ia(l_dof)+1
                               else 
                                  if ( m_dof >= l_dof ) then
                                     ic = aux_ia(l_dof)
                                     dof_graph%ja(ic) = m_dof
                                     aux_ia(l_dof) = aux_ia(l_dof)+1
                                  end if
                               end if
                            end if
                         end do
                      end do
                   end if
                end do
                !end do
                if (.not.fe_space%static_condensation) then  ! interface-interior
                   iprob = fe_space%finite_elements(jelem)%problem
                   nvapb = dof_descriptor%prob_block(jblock,iprob)%nd1
                   do idof = fe_space%vef2dof(iblock)%p(iobje), fe_space%vef2dof(iblock)%p(iobje+1)-1
                      l_dof = fe_space%vef2dof(iblock)%l(idof,1)
                      l_var = fe_space%vef2dof(iblock)%l(idof,2)
                      l_mat = fe_space%vef2dof(iblock)%l(idof,3)
                      do ivars = 1, nvapb
                         k_var = dof_descriptor%prob_block(jblock,iprob)%a(ivars)
                         m_var = dof_descriptor%problems(iprob)%p%l2g_var(k_var)
                         m_mat = fe_space%finite_elements(jelem)%continuity(m_var)
                         if ( dof_descriptor%dof_coupl(l_var, m_var) == 1 .and. l_mat == m_mat ) then                
                            if ( .not. dof_graph%symmetric_storage ) then
                               do inode = fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(jobje), &
                                    & fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(jobje+1)-1
                                  l_node = fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%l(inode)
                                  m_dof = fe_space%finite_elements(jelem)%elem2dof(l_node,k_var)
                                     ic = aux_ia(l_dof)
                                     dof_graph%ja(ic) = m_dof
                                     aux_ia(l_dof) = aux_ia(l_dof)+1
                               end do
                            else 
                               do inode = fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(jobje), &
                                    & fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(jobje+1)-1
                                  l_node = fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%l(inode)
                                  m_dof = fe_space%finite_elements(jelem)%elem2dof(l_node,k_var)
                                  if ( m_dof >= l_dof ) then
                                     ic = aux_ia(l_dof)
                                     dof_graph%ja(ic) = m_dof
                                     aux_ia(l_dof) = aux_ia(l_dof)+1
                                  end if
                               end do
                            end if
                         end if
                      end do
                   end do
                end if
             end if
          end do
          call visited%free
       end if
    end do

  end subroutine list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity

  !*********************************************************************************
  ! Count NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
  ! both interior and interface nodes.
  !*********************************************************************************
  subroutine count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dof_descriptor, trian, fe_space, dof_graph )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_descriptor_t), intent(in)               :: dof_descriptor
    type(triangulation_t), intent(in)         :: trian 
    type(fe_space_t), intent(in)                 :: fe_space 
    type(graph_t), intent(inout)                :: dof_graph

    ! Local variables
    integer(ip) :: g_var, ielem, inode, int_i, iobje, iprob, ivars, jdof, jnode, job_g
    integer(ip) :: jobje, jvars, k_var, l_dof, l_mat, l_node, l_var, m_dof, m_mat
    integer(ip) :: m_node, m_var, nvapbi, nvapbj

    ! As commented for elem2dof, static condensation is false for dG, by construction of the 
    ! fem space.
    if (.not.fe_space%static_condensation) then
       do ielem  = 1, trian%num_elems
          iobje = trian%elems(ielem)%num_vefs+1
          iprob = fe_space%finite_elements(ielem)%problem
          nvapbi = dof_descriptor%prob_block(iblock,iprob)%nd1 
          do ivars = 1, nvapbi
             l_var = dof_descriptor%prob_block(iblock,iprob)%a(ivars)
             g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
             ! Interior - interior 
             nvapbj = dof_descriptor%prob_block(jblock,iprob)%nd1 
             do jvars = 1, nvapbj
                k_var = dof_descriptor%prob_block(jblock,iprob)%a(jvars)
                m_var = dof_descriptor%problems(iprob)%p%l2g_var(k_var)
                if ( dof_descriptor%dof_coupl(g_var,m_var) == 1 ) then
                   do inode = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje), &
                        & fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje+1)-1
                      l_node = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%l(inode)
                      l_dof = fe_space%finite_elements(ielem)%elem2dof(l_node,l_var)
                      if ( .not. dof_graph%symmetric_storage ) then
                         dof_graph%ia(l_dof+1) =  dof_graph%ia(l_dof+1) &
                              &  + fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%p(iobje+1) &
                              & - fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%p(iobje)
                      else  
                         do jnode = fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%p(iobje), &
                              & fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%p(iobje+1)-1
                            m_node = fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%l(jnode)
                            m_dof = fe_space%finite_elements(ielem)%elem2dof(m_node,k_var)
                            if ( m_dof >= l_dof ) then
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1
                            end if
                         end do
                      end if
                   end do
                end if
             end do
             l_mat = fe_space%finite_elements(ielem)%continuity(g_var)
             if ( l_mat /= 0 ) then
                ! Interior - interface 
                do jobje = 1, trian%elems(ielem)%num_vefs
                   job_g = trian%elems(ielem)%vefs(jobje)
                   do jdof = fe_space%vef2dof(jblock)%p(job_g), fe_space%vef2dof(jblock)%p(job_g+1)-1
                      m_dof = fe_space%vef2dof(jblock)%l(jdof,1)
                      m_var = fe_space%vef2dof(jblock)%l(jdof,2)   
                      m_mat = fe_space%vef2dof(jblock)%l(jdof,3)                      
                      if ( dof_descriptor%dof_coupl(g_var,m_var) == 1 .and. l_mat == m_mat ) then
                         do inode = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje), &
                              & fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje+1)-1
                            l_node = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%l(inode)
                            l_dof = fe_space%finite_elements(ielem)%elem2dof(l_node,l_var)
                            if ( .not. dof_graph%symmetric_storage ) then
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1 
                            else if ( m_dof >= l_dof ) then
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end do
       end do
    end if

  end subroutine count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity

  !*********************************************************************************
  ! List NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
  ! both interior and interface nodes.
  !*********************************************************************************
  subroutine list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dof_descriptor, trian, fe_space, dof_graph, aux_ia )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_descriptor_t), intent(in)               :: dof_descriptor
    type(triangulation_t), intent(in)         :: trian 
    type(fe_space_t), intent(in)                 :: fe_space 
    type(graph_t), intent(inout)              :: dof_graph
    integer(ip), intent(inout)                  :: aux_ia(:) 

    ! Local variables
    integer(ip) :: g_var, ielem, inode, iobje, iprob, ivars, jdof, jnode, job_g
    integer(ip) :: jobje, jvars, k_var, l_dof, l_mat, l_node, l_var, m_dof, m_mat
    integer(ip) :: m_node, m_var, nvapbi, nvapbj, i, ic

    if (.not.fe_space%static_condensation) then   
       do ielem  = 1, trian%num_elems
          iobje = trian%elems(ielem)%num_vefs+1
          iprob = fe_space%finite_elements(ielem)%problem
          nvapbi = dof_descriptor%prob_block(iblock,iprob)%nd1  
          do ivars = 1, nvapbi
             !l_var = g2l(ivars,iprob)
             l_var = dof_descriptor%prob_block(iblock,iprob)%a(ivars)
             g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
             ! Interior - interior (inside element)
             nvapbj = dof_descriptor%prob_block(jblock,iprob)%nd1  
             do jvars = 1, nvapbj
                k_var = dof_descriptor%prob_block(jblock,iprob)%a(jvars)
                m_var = dof_descriptor%problems(iprob)%p%l2g_var(k_var)
                if ( dof_descriptor%dof_coupl(g_var,m_var) == 1 ) then
                   do inode = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje), &
                        & fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje+1)-1
                      l_node = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%l(inode)
                      l_dof = fe_space%finite_elements(ielem)%elem2dof(l_node,l_var)
                      if ( .not. dof_graph%symmetric_storage ) then
                         do jnode = fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%p(iobje), &
                              & fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%p(iobje+1)-1
                            m_node = fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%l(jnode)
                            m_dof = fe_space%finite_elements(ielem)%elem2dof(m_node,k_var)
                            i= aux_ia(l_dof)
                            dof_graph%ja(i) = m_dof
                            aux_ia(l_dof) = aux_ia(l_dof)+1
                         end do
                      else ! ltype == csr_symm 
                         do jnode = fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%p(iobje), &
                              & fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%p(iobje+1)-1
                            m_node = fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%l(jnode)
                            m_dof = fe_space%finite_elements(ielem)%elem2dof(m_node,k_var)
                            if ( m_dof >= l_dof ) then
                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof)+1
                            end if
                         end do
                      end if
                   end do
                end if
             end do
             if ( fe_space%finite_elements(ielem)%continuity(g_var) /= 0 ) then
                ! Interior - border (inside element)
                do jobje = 1, trian%elems(ielem)%num_vefs
                   job_g = trian%elems(ielem)%vefs(jobje)
                   do jdof = fe_space%vef2dof(jblock)%p(job_g), fe_space%vef2dof(jblock)%p(job_g+1)-1
                      m_dof = fe_space%vef2dof(jblock)%l(jdof,1)
                      m_var = fe_space%vef2dof(jblock)%l(jdof,2)                         
                      if ( dof_descriptor%dof_coupl(g_var,m_var) == 1 ) then
                         do inode = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje), &
                              & fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje+1)-1
                            l_node = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%l(inode)
                            l_dof = fe_space%finite_elements(ielem)%elem2dof(l_node,l_var)
                            if ( .not. dof_graph%symmetric_storage ) then
                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof)+1
                            else if ( m_dof >= l_dof ) then
                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof)+1
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end do
       end do
    end if
  end subroutine list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity

  !*********************************************************************************
  ! Count NNZ (number of nonzero entries) for DOFs on the element interior (for dG) being 
  ! coupled due to integration on faces. *** This part requires more ellaboration for cdG
  ! generalization***
  !*********************************************************************************
  ! Note: Here we take a face in which we want to integrate dG terms for pairs of unknowns
  ! and couple all DOFs a la DG. Note that we are not using AT ALL the continuity value,
  ! since the integration in faces is driven by the faces selected for integration, which
  ! are being built accordingly to what one wants to do (cG, dG, dG for jump of cont value, etc.).
  ! One could think that it could happen that two elements K1 and K2 with two different 
  ! values of continuity that share a face where we want to integrate could put more than
  ! once the coupling among two nodes. It can never happen AS SOON AS one never creates
  ! an integration face between two elements with same continuity value (not 0), which is the
  ! expected usage. 
  ! *** We could put an assert about it when creating the integration list.
  !*********************************************************************************
  subroutine count_nnz_all_dofs_vs_all_dofs_by_face_integration ( iblock, jblock, dof_descriptor, trian, fe_space, dof_graph )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_descriptor_t), intent(in)               :: dof_descriptor
    type(triangulation_t), intent(in)         :: trian 
    type(fe_space_t), intent(in)                 :: fe_space 
    type(graph_t), intent(inout)              :: dof_graph

    ! Local variables
    integer(ip) :: count, g_var, i, ielem, iface, inode, iobje, iprob, ivars, jelem
    integer(ip) :: jnode, jprob, jvars, k_var, knode, l_dof, l_faci, l_facj, l_node
    integer(ip) :: l_var, m_dof, m_var, m_node, nnode, nvapbi, nvapbj

    ! Loop over all interior faces (boundary faces do not include additional coupling)
    do iface = 1,fe_space%num_interior_faces
       iobje = fe_space%interior_faces(iface)%face_vef
       assert ( trian%vefs(iobje)%num_elems_around == 2 ) 
       do i=1,2
          ielem = trian%vefs(iobje)%elems_around(i)
          jelem = trian%vefs(iobje)%elems_around(3-i)
          l_faci = local_position(fe_space%interior_faces(iface)%face_vef,trian%elems(ielem)%vefs, &
               & trian%elems(ielem)%num_vefs )
          l_facj = local_position(fe_space%interior_faces(iface)%face_vef,trian%elems(jelem)%vefs, &
               & trian%elems(jelem)%num_vefs )
          iprob = fe_space%finite_elements(ielem)%problem
          jprob = fe_space%finite_elements(jelem)%problem
          nvapbi = dof_descriptor%prob_block(iblock,iprob)%nd1 
          nvapbj = dof_descriptor%prob_block(iblock,jprob)%nd1 
          do ivars = 1, nvapbi
             l_var = dof_descriptor%prob_block(iblock,iprob)%a(ivars)
             g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
             do jvars = 1, nvapbj
                k_var = dof_descriptor%prob_block(iblock,jprob)%a(jvars)
                m_var = dof_descriptor%problems(jprob)%p%l2g_var(k_var)
                if ( dof_descriptor%dof_coupl(g_var,m_var) == 1 ) then
                   if ( .not. dof_graph%symmetric_storage ) then
                      ! Couple all DOFs in ielem with face DOFs in jelem and viceversa (i=1,2)
                      nnode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1) &
                           &  -fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj)
                      do inode = 1, fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%nnode
                         l_dof = fe_space%finite_elements(ielem)%elem2dof(inode,l_var)
                         if ( l_dof /= 0 ) then
                            dof_graph%ia(l_dof+1) = dof_graph%ia(l_dof+1) &
                                 & + nnode
                         end if
                      end do
                      nnode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%nnode - nnode
                      assert ( nnode > 0)
                      ! Couple all face DOFs in ielem with DOFs in jelem NOT in the face and viceversa (i=1,2)
                      do inode = fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%p(l_faci), &
                           &     fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%p(l_faci+1)-1
                         l_node = fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%l(inode)
                         l_dof = fe_space%finite_elements(ielem)%elem2dof(l_node,l_var)
                         if ( l_dof /= 0 ) then
                            dof_graph%ia(l_dof+1) = dof_graph%ia(l_dof+1) &
                                 & + nnode
                         end if
                      end do
                   else ! ltype == csr_symm 
                      ! Couple all DOFs in ielem with face DOFs in jelem and viceversa (i=1,2)
                      do inode = 1, fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%nnode
                         l_dof = fe_space%finite_elements(ielem)%elem2dof(inode,l_var)
                         if ( l_dof /= 0 ) then
                            do jnode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj), &
                                 & fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1
                               m_node = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(jnode)
                               m_dof = fe_space%finite_elements(jelem)%elem2dof(m_node,k_var)
                               if ( l_dof /= 0 .and. m_dof >= l_dof ) then
                                  dof_graph%ia(l_dof+1) = &
                                       & dof_graph%ia(l_dof+1) + 1
                               end if
                            end do
                         end if
                      end do
                      ! Couple all face DOFs in ielem with DOFs in jelem NOT in the face and viceversa (i=1,2)
                      count = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(count)
                      end if
                      do jnode = 1, fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%nnode
                         if ( jnode == knode) then
                            count = count+1
                            if (count <= fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(count)
                            end if
                         else
                            m_node = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(jnode)
                            m_dof = fe_space%finite_elements(jelem)%elem2dof(m_node,k_var)

                            do inode = fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%p(l_faci), &
                                 &     fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%p(l_faci+1)-1
                               l_node = fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%l(inode)
                               l_dof = fe_space%finite_elements(ielem)%elem2dof(l_node,l_var)
                               if ( m_dof >= l_dof ) then
                                  !write (*,*) 'INSERTION DUE TO COUPLING BY FACE(ielem)-INTERIOR(jelem)'
                                  !write (*,*) 'ielem,jelem,iobje,l_faci,l_facj',ielem,jelem,iobje,l_faci,l_facj
                                  !write (*,*) 'ielem,jelem'
                                  !write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
                                  !write (*,*) 'NOW',l_dof,'COUPLED TO',dof_graph%ia(l_dof+1) + 1
                                  !write (*,*) 'CHECK COLISION:count,knode,jnode',count,knode,jnode
                                  dof_graph%ia(l_dof+1) = &
                                       & dof_graph%ia(l_dof+1) + 1
                               end if
                            end do
                         end if
                      end do

                   end if
                end if
             end do
          end do
       end do
    end do
  end subroutine count_nnz_all_dofs_vs_all_dofs_by_face_integration

  !*********************************************************************************
  ! Count NNZ (number of nonzero entries) for DOFs on the element interior (for dG) being 
  ! coupled due to integration on faces. *** This part requires more ellaboration for cdG
  ! generalization***
  !*********************************************************************************
  subroutine list_nnz_all_dofs_vs_all_dofs_by_face_integration ( iblock, jblock, dof_descriptor, trian, fe_space, dof_graph, aux_ia )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_descriptor_t), intent(in)               :: dof_descriptor
    type(triangulation_t), intent(in)         :: trian 
    type(fe_space_t), intent(in)                 :: fe_space 
    type(graph_t), intent(inout)              :: dof_graph
    integer(ip), intent(inout)                  :: aux_ia(:) 

    ! Local variables
    integer(ip) :: count, g_var, i, ielem, iface, inode, iobje, iprob, ivars, jelem
    integer(ip) :: jnode, jprob, jvars, k_var, knode, l_dof, l_faci, l_facj, l_node
    integer(ip) :: l_var, m_dof, m_var, m_node, nnode, nvapbi, nvapbj, ic

    ! Loop over all interior faces (boundary faces do not include additional coupling)
    do iface = 1, fe_space%num_interior_faces
       iobje = fe_space%interior_faces(iface)%face_vef
       assert ( trian%vefs(iobje)%num_elems_around == 2 ) 
       do i=1,2
          ielem = trian%vefs(iobje)%elems_around(i)
          jelem = trian%vefs(iobje)%elems_around(3-i)
          l_faci = local_position(iobje,trian%elems(ielem)%vefs, &
               & trian%elems(ielem)%num_vefs )
          l_facj = local_position(iobje,trian%elems(jelem)%vefs, &
               & trian%elems(jelem)%num_vefs )
          iprob = fe_space%finite_elements(ielem)%problem
          jprob = fe_space%finite_elements(jelem)%problem
          nvapbi = dof_descriptor%prob_block(iblock,iprob)%nd1 
          nvapbj = dof_descriptor%prob_block(iblock,jprob)%nd1 
          do ivars = 1, nvapbi
             l_var = dof_descriptor%prob_block(iblock,iprob)%a(ivars)
             g_var = dof_descriptor%problems(iprob)%p%l2g_var(l_var)
             do jvars = 1, nvapbj
                k_var = dof_descriptor%prob_block(iblock,jprob)%a(jvars)
                m_var = dof_descriptor%problems(jprob)%p%l2g_var(k_var)
                if ( dof_descriptor%dof_coupl(g_var,m_var) == 1 ) then
                   if ( .not. dof_graph%symmetric_storage ) then
                      ! Couple all DOFs in ielem with face DOFs in jelem and viceversa (i=1,2)
                      do inode = 1, fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%nnode
                         l_dof = fe_space%finite_elements(ielem)%elem2dof(inode,l_var)
                         if ( l_dof /= 0 ) then 
                            !                         do jnode = fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(l_facj), &
                            !                              &     fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%p(l_facj+1)-1 
                            !                            m_node = fe_space%finite_elements(ielem)%nodes_per_vef(k_var)%p%l(jnode)
                            do jnode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj), &
                                 & fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1
                               m_node = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(jnode)
                               !                            m_node = fe_space%finite_elements(jelem)%nodes_per_vef(k_var)%p%l(jnode)
                               m_dof = fe_space%finite_elements(jelem)%elem2dof(m_node,k_var)
                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof) + 1
                            end do
                         end if
                      end do
                      ! Couple all face DOFs in ielem with DOFs in jelem NOT in the face and viceversa (i=1,2)
                      count = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(count)
                      end if
                      do jnode = 1, fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%nnode
                         if ( jnode == knode) then
                            count = count+1
                            if (count <= fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(count)
                            end if
                         else
                            m_node = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(jnode)
                            m_dof = fe_space%finite_elements(jelem)%elem2dof(m_node,k_var)
                            do inode = fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%p(l_faci), &
                                 &     fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%p(l_faci+1)-1
                               l_node = fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%l(inode)
                               l_dof = fe_space%finite_elements(ielem)%elem2dof(l_node,l_var)
                               if (l_dof /= 0 ) then
                                  ic = aux_ia(l_dof)
                                  dof_graph%ja(ic) = m_dof
                                  aux_ia(l_dof) = aux_ia(l_dof) + 1
                               end if
                            end do
                         end if
                      end do
                   else ! ltype == csr_symm 
                      ! Couple all DOFs in ielem with face DOFs in jelem and viceversa (i=1,2)
                      do inode = 1, fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%nnode
                         l_dof = fe_space%finite_elements(ielem)%elem2dof(inode,l_var)
                         if (l_dof /= 0 ) then
                            do jnode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj), &
                                 & fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1
                               m_node = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(jnode)
                               m_dof = fe_space%finite_elements(jelem)%elem2dof(m_node,k_var)
                               if ( m_dof >= l_dof ) then
                                  !write (*,*) 'VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV'
                                  !write (*,*) 'INSERTION DUE TO COUPLING BY (ielem)-FACE(jelem)'
                                  !write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
                                  !write (*,*) 'IN POSITION',aux_ia(l_dof)
                                  !write (*,*) 'VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV'
                                  ic = aux_ia(l_dof)
                                  dof_graph%ja(ic) = m_dof
                                  aux_ia(l_dof) = aux_ia(l_dof) + 1
                               end if
                            end do
                         END if
                      end do
                      ! Couple all face DOFs in ielem with DOFs in jelem NOT in the face and viceversa (i=1,2)
                      count = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(count)
                      end if
                      do jnode = 1, fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%nnode
                         if ( jnode == knode) then
                            count = count+1
                            if (count <= fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(count)
                            end if
                         else
                            m_node = fe_space%finite_elements(jelem)%reference_element_vars(k_var)%p%ntxob%l(jnode)
                            m_dof = fe_space%finite_elements(jelem)%elem2dof(m_node,k_var)
                            do inode = fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%p(l_faci), &
                                 &     fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%p(l_faci+1)-1
                               l_node = fe_space%finite_elements(ielem)%reference_element_vars(l_var)%p%ntxob%l(inode)
                               l_dof = fe_space%finite_elements(ielem)%elem2dof(l_node,l_var)
                               if ( l_dof /= 0 .and. m_dof >= l_dof ) then
                                  !write (*,*) 'KKKKKKKKKKKKKK'
                                  !write (*,*) 'INSERTION DUE TO COUPLING BY FACE(ielem)-INTERIOR(jelem)'
                                  !write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
                                  !write (*,*) 'IN POSITION',aux_ia(l_dof)
                                  !write (*,*) 'KKKKKKKKKKKKKK'
                                  ic = aux_ia(l_dof)
                                  dof_graph%ja(ic) = m_dof
                                  aux_ia(l_dof) = aux_ia(l_dof) + 1
                               end if
                            end do
                         end if
                      end do
                   end if
                end if
             end do
          end do
       end do
    end do

  end subroutine list_nnz_all_dofs_vs_all_dofs_by_face_integration

  !*********************************************************************************
  ! Auxiliary function that returns the position of key in list
  !*********************************************************************************
  integer(ip) function local_position(key,list,size)
    implicit none
    integer(ip) :: key, size, list(size)

    do local_position = 1,size
       if ( list(local_position) == key) exit
    end do
    assert ( local_position < size + 1 )

  end function local_position

end module fe_space_names

