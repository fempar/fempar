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
# include "debug.i90"
module fem_space_names
  ! Modules
  use types
  use memor
  use array_names
  use fem_mesh_names
  use hash_table_names
!  use elmat_names
!  use elvec_names
  use integration_names
  use fem_space_types
#ifdef memcheck
  use iso_c_binding
#endif
  implicit none
  private

  ! Information of each element of the FE space
  type fem_element

     integer(ip)                   :: prob          ! Problem to be solved
     integer(ip)                   :: nint          ! Number of different interpolations per element
     integer(ip)                   :: ldof          ! Dofs per element

     integer(ip)     , allocatable :: elem2dof(:,:) ! Map from elem to dof
     integer(ip)     , allocatable :: iv(:)         ! Interpolation of each variable
     integer(ip)     , allocatable :: p_nod(:)      ! Pointer to first DOF on elemental matrix/vector
     integer(ip)     , allocatable :: jvars(:,:)    ! Variable identifier for each node

     real(rp)        , allocatable :: unkno(:,:,:)  ! Values of the solution on the nodes of the elem    
   
     type(array_rp1) , allocatable :: bc_value(:)   ! Boundary Condition values
     type(fem_fixed_info), pointer :: p_geo_info => NULL() ! Interpolation info of the geometry
     type(fem_fixed_info_pointer), allocatable :: f_inf(:) ! Interpolation info of the FE space
    
     type(array_rp2), pointer :: p_mat ! Pointer to the elemental matrix
     type(array_rp1), pointer :: p_vec ! Pointer to the elemental vector

     type(vol_integ_pointer)     , allocatable :: integ(:)  ! Pointer to integration parameters
  end type fem_element
  
  ! Information relative to the faces
  type fem_face
     integer(ip)               :: nei_elem(2)       ! Neighbor elements
     integer(ip)               :: pos_elem(2)       ! Face pos in element
     integer(ip)               :: subface(2)        ! It can be a portion of one of the elems face
     integer(ip)               :: refinement_level(2)

     type(array_rp2), pointer  :: p_mat ! Pointer to the elemental matrix

     !type(array_rp2), pointer  :: aux_mat(2)        ! Pointer to integration face matrices
     type(array_rp1), pointer  :: p_vec   ! Pointer to face integration vector
     type(face_integ_pointer), allocatable :: integ(:)  ! Pointer to face integration

     !integer(ip) , allocatable :: o2n(:)            ! permutation of the gauss points in elem2
     type(array_ip1), allocatable:: o2n(:)           ! permutation of the gauss points in elem2
  end type fem_face

  ! Global information of the fem_space
  type fem_space  

     !integer(ip)                    :: nelem             ! Number of elements
     !integer(ip)                    :: dim_space         ! Dimension of the space
     !type(dof_handler)              :: f_dof_handler
     integer(ip)                        :: kfl_scond        ! Flag for static condensation 
                                                            !(0=deactivated,1=activated)
     ! Elements and faces
     type(fem_element)    , allocatable :: lelem(:)         ! List of FEs
     type(fem_face)       , allocatable :: lface(:)         ! List of active faces

     type(fem_mesh)           , pointer :: o_mesh => NULL() ! Mesh of objects
     type(fem_mesh)           , pointer :: g_mesh => NULL() ! Mesh of geometrical nodes (as now)

     ! Elemental matrices and vector
     type (hash_table_ip_ip)            :: ht_pos_elmat
     type (hash_table_ip_ip)            :: ht_pos_elvec
     type(array_rp2), allocatable       :: lelmat(:)
     type(array_rp1), allocatable       :: lelvec(:)

     integer(ip)                        :: cur_elmat
     integer(ip)                        :: cur_elvec

     ! Integrator
     type (hash_table_ip_ip)            :: ht_pos_vol_integ
     type (hash_table_ip_ip)            :: ht_pos_face_integ
     type(vol_integ)      , allocatable :: lvoli(:)        
     type(face_integ)     , allocatable :: lfaci(:)
     integer(ip)                        :: cur_lvoli
     integer(ip)                        :: cur_lfaci

     ! Element common information
     type (hash_table_ip_ip)            :: ht_elem_info
     type (fem_fixed_info), allocatable :: lelem_info(:)
     integer(ip)                         :: cur_elinf

  end type fem_space

  interface fem_space_fe_list_create
     module procedure fem_space_fe_list_create_one_int, fem_space_fe_list_create_var_int
  end interface fem_space_fe_list_create

  ! Types
  public :: fem_space, fem_element, fem_face,                      &
       &     memalloc, memrealloc, memfreep, memmovealloc

  ! Functions
  public :: fem_space_create, fem_space_fe_list_create, fem_element_print, fem_space_free,          &
            get_p_faces, fem_space_fe_list_create_one_int_p
  
!!$   integer(ip) :: P3_connec(3,14) = reshape((/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
!!$        & 1, 2, 0, 2, 3, 0, 3, 1, 0, 1, 4, 0, 2, 4, 0, 3, 4, 0, &
!!$        & 1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4  /),(/3,14/))
!!$   
!!$   integer(ip) :: Q3_connec(4,26) = reshape((/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
!!$        & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
!!$        & 1, 2, 0, 0, 2, 3, 0, 0, 3, 4, 0, 0, 4, 1, 0, 0, 1, 5, 0, 0, 2, 6, 0, 0, &
!!$        & 3, 7, 0, 0, 4, 8, 0, 0, 5, 6, 0, 0, 6, 7, 0, 0, 7, 8, 0, 0, 8, 5, 0, 0, &
!!$        & 1, 2, 4, 3, 1, 2, 5, 6, 1, 5, 4, 8, 2, 6, 3, 7, 3, 4, 7, 8, 5, 6, 8, 7 /),(/4,26/))
   
contains

  !==================================================================================================
  ! Allocation of variables in fem_space according to the values in o_mesh and g_mesh
  subroutine fem_space_create( fspac, o_mesh, g_mesh)
    implicit none
    type(fem_space)  , intent(inout)        :: fspac
    type(fem_mesh)   , intent(in), target   :: o_mesh, g_mesh

    integer(ip) :: istat

    allocate(fspac%lelem( g_mesh%nelem), stat=istat)
    if(istat/=0) then 
       write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
       write(6,'(a,i20)') 'Error code: ', istat
       write(6,'(a)') 'Caller routine: fem_space::fem_space_create::lelem'
    end if

    !  Initialization of pointers to topological and geometrical mesh
    fspac%o_mesh => o_mesh
    fspac%g_mesh => g_mesh

    ! Allocation of elemental matrix and vector parameters
    call fspac%ht_pos_elmat%init(ht_length)
    fspac%cur_elmat = 1
    allocate(fspac%lelmat(max_elmat), stat=istat)
    if(istat/=0) then 
       write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
       write(6,'(a,i20)') 'Error code: ', istat
       write(6,'(a)') 'Caller routine: fem_space::fem_space_create::lelmat'
    end if
    call fspac%ht_pos_elvec%init(ht_length)
    fspac%cur_elvec = 1
    allocate(fspac%lelvec(max_elmat), stat=istat)
    if(istat/=0) then 
       write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
       write(6,'(a,i20)') 'Error code: ', istat
       write(6,'(a)') 'Caller routine: fem_space::fem_space_create::lelvec'
    end if

    ! Initialization of volume and face integrators parameters
    call fspac%ht_pos_vol_integ%init(ht_length)
    call fspac%ht_pos_face_integ%init(ht_length)
    fspac%cur_lvoli = 1
    fspac%cur_lfaci = 1
    allocate(fspac%lvoli(max_elmat), stat=istat)
    if(istat/=0) then 
       write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
       write(6,'(a,i20)') 'Error code: ', istat
       write(6,'(a)') 'Caller routine: fem_space::fem_space_create::lvoli'
    end if

    ! Initialization of element fixed info parameters
    call fspac%ht_elem_info%init(ht_length)
    fspac%cur_elinf = 1
    allocate(fspac%lelem_info(max_elinf), stat=istat)
    if(istat/=0) then 
       write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
       write(6,'(a,i20)') 'Error code: ', istat
       write(6,'(a)') 'Caller routine: fem_space::fem_space_create::lelem_info'
    end if
  end subroutine fem_space_create

  !==================================================================================================
  ! Fill the fem_space assuming that all elements are of type f_type and order f_order
  subroutine fem_space_fe_list_create_one_int( fspac, f_type, f_order, dim, tdim, prob_code,        &
       prob_nunk, prob_list, prob_list_nunk, khie, scond  )
    implicit none
    ! Parameters
    type(fem_space), intent(inout),target  :: fspac
    integer(ip)          , intent(in) :: f_type, f_order, dim, tdim
    integer(ip), optional, intent(in) :: prob_code, prob_nunk
    integer(ip), optional, intent(in) :: prob_list(:), prob_list_nunk(:)
    integer(ip), optional, intent(in) :: khie
    integer(ip), optional, intent(in) :: scond

    ! Local variables
    integer(ip) :: nunk, v_key, ltype(2), khir, nnode, nunk_tot
    integer(ip) :: ielem, istat
    integer(ip) :: pos_elinf, pos_elmat, pos_elvec, pos_voint
    logical(ip) :: created
    
    ! Hierarchical flag
    if (present(khie)) then
       khir = khie
    else
       khir = 0
    end if
    
    ! Static condensation flag
    if (present(scond)) then
       fspac%kfl_scond = scond
    else
       fspac%kfl_scond = scond_off
    end if

    ! nunk_tot = total amount of unknowns
    if ( present(prob_list_nunk) ) then
       nunk_tot=prob_list_nunk(1)
       do ielem=2, fspac%g_mesh%nelem
          nunk_tot=max(nunk_tot,prob_list_nunk(ielem))
       end do
    else
       nunk_tot = prob_nunk
    end if

    ! Loop over elements
    do ielem = 1, fspac%g_mesh%nelem
       ! Assign type of problem to ielem
       if ( present(prob_code) .and. present(prob_nunk) ) then
          fspac%lelem(ielem)%prob = prob_code
          nunk = prob_nunk
       else if ( present(prob_list) .and. present(prob_list_nunk) ) then
          fspac%lelem(ielem)%prob = prob_list(ielem)
          nunk = prob_list_nunk(ielem)
       else
          write(6,*) 'Error in fem_space_fe_list_create'
          assert(0==1)
       end if

       ! One only interpolation per element
       fspac%lelem(ielem)%nint = 1

       ! All the variables have the same interpolation
       call memalloc(nunk_tot,fspac%lelem(ielem)%iv,__FILE__,__LINE__)
       fspac%lelem(ielem)%iv(:) = 1

       ! Allocate the fixed info array
       call memalloc(1,fspac%lelem(ielem)%f_inf,__FILE__,__LINE__)

       ! Assign pointer to interpolation fixed information
       v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
       call fspac%ht_elem_info%put(key=v_key,val=fspac%cur_elinf,stat=istat)
       if ( istat == now_stored) then 
          ! Create fixed info if not constructed
          call fem_element_fixed_info_create(fspac%lelem_info(fspac%cur_elinf),f_type,              &
               &                             f_order,dim,created)
          assert(created)
          pos_elinf = fspac%cur_elinf
          fspac%cur_elinf = fspac%cur_elinf + 1
       else if ( istat == was_stored ) then
          call fspac%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          assert ( istat == key_found )
       end if
       fspac%lelem(ielem)%f_inf(1)%p => fspac%lelem_info(pos_elinf)

       ! Assign pointer to geometrical fixed information (assumed to be of order 1)
       v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       call fspac%ht_elem_info%put(key=v_key,val=fspac%cur_elinf,stat=istat)
       if ( istat == now_stored) then 
          ! Create fixed info if not constructed
          call fem_element_fixed_info_create(fspac%lelem_info(fspac%cur_elinf),f_type,              &
               &                             1,dim,created)
          assert(created)
          pos_elinf = fspac%cur_elinf
          fspac%cur_elinf = fspac%cur_elinf + 1
       else if ( istat == was_stored ) then
          call fspac%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          assert ( istat == key_found )
       end if
       fspac%lelem(ielem)%p_geo_info => fspac%lelem_info(pos_elinf)

       ! Compute number of DOFs in the elemental matrix associated to ielem
       if (fspac%kfl_scond == scond_on) then
          nnode = fspac%lelem(ielem)%f_inf(1)%p%nnode-fspac%lelem(ielem)%f_inf(1)%p%nodes_obj(dim+1)
       else
          nnode = fspac%lelem(ielem)%f_inf(1)%p%nnode
       end if
       fspac%lelem(ielem)%ldof = nunk*nnode

       ! Assign pointer to p_mat and p_vec in ielem
       call fspac%ht_pos_elmat%put(key=nunk*nnode,val=fspac%cur_elmat,stat=istat)
       call fspac%ht_pos_elvec%put(key=nunk*nnode,val=fspac%cur_elvec,stat=istat)
       if ( istat == now_stored ) then
          call elmat_alloc ( 1, 1, nunk*nnode,fspac%lelmat(fspac%cur_elmat), scal )
          pos_elmat = fspac%cur_elmat
          fspac%cur_elmat = fspac%cur_elmat + 1
          call elvec_alloc (1,nunk*nnode,fspac%lelvec(fspac%cur_elvec), scal )
          pos_elvec = fspac%cur_elvec
          fspac%cur_elvec = fspac%cur_elvec + 1
       else if ( istat == was_stored ) then
          call fspac%ht_pos_elmat%get(key=nunk*nnode,val=pos_elmat,stat=istat)
          assert ( istat == key_found )
          call fspac%ht_pos_elvec%get(key=nunk*nnode,val=pos_elvec,stat=istat)
          assert ( istat == key_found )
       end if       
       fspac%lelem(ielem)%p_mat => fspac%lelmat(pos_elmat)
       fspac%lelem(ielem)%p_vec => fspac%lelvec(pos_elvec)
       
       ! Allocate elem2dof, unkno, p_nod, integ and jvars
       call memalloc( nnode, nunk_tot, fspac%lelem(ielem)%elem2dof, __FILE__,__LINE__ )
       fspac%lelem(ielem)%elem2dof = 0
       call memalloc( nnode, nunk_tot, tdim,fspac%lelem(ielem)%unkno, __FILE__,__LINE__)
       fspac%lelem(ielem)%unkno = 0.0_rp
       call memalloc( nnode+1, fspac%lelem(ielem)%p_nod, __FILE__,__LINE__)
       fspac%lelem(ielem)%p_nod = 0
       call memalloc(1,fspac%lelem(ielem)%integ,__FILE__,__LINE__)
       call memalloc( nnode, nunk_tot, fspac%lelem(ielem)%jvars, __FILE__,__LINE__ )
       fspac%lelem(ielem)%jvars = 0

       ! Assign pointers to volume integration
       ltype(1) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order
       ltype(2) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       v_key    = (max_ndime+1)*(max_FE_types+1)*(max_order) * ltype(1) + ltype(2)
       call fspac%ht_pos_vol_integ%put(key=v_key, val=fspac%cur_lvoli, stat = istat)
       if ( istat == now_stored ) then
          call integ_create(f_type,f_type,dim,1,f_order,fspac%lvoli(fspac%cur_lvoli),khie = khir)
          pos_voint       = fspac%cur_lvoli
          fspac%cur_lvoli = fspac%cur_lvoli + 1
       else if ( istat == was_stored ) then
          call fspac%ht_pos_vol_integ%get(key=v_key,val=pos_voint,stat=istat)
          assert ( istat == key_found )
       end if
       fspac%lelem(ielem)%integ(1)%p => fspac%lvoli(pos_voint)
       
    end do

  end subroutine fem_space_fe_list_create_one_int

  !==================================================================================================
  ! Fill the fem_space assuming that all elements are of type f_type and order f_order
  subroutine fem_space_fe_list_create_one_int_p( fspac, f_type, f_order, dim, tdim, prob_code,        &
       prob_nunk, prob_list, prob_list_nunk, khie, scond  )
    implicit none
    ! Parameters
    type(fem_space), intent(inout),target  :: fspac
    integer(ip)          , intent(in) :: f_type, f_order(:), dim, tdim
    integer(ip), optional, intent(in) :: prob_code, prob_nunk
    integer(ip), optional, intent(in) :: prob_list(:), prob_list_nunk(:)
    integer(ip), optional, intent(in) :: khie
    integer(ip), optional, intent(in) :: scond

    ! Local variables
    integer(ip) :: nunk, v_key, ltype(2), khir, nnode, nunk_tot
    integer(ip) :: ielem, istat
    integer(ip) :: pos_elinf, pos_elmat, pos_elvec, pos_voint
    logical(ip) :: created
    
    ! Hierarchical flag
    if (present(khie)) then
       khir = khie
    else
       khir = 0
    end if
    
    ! Static condensation flag
    if (present(scond)) then
       fspac%kfl_scond = scond
    else
       fspac%kfl_scond = scond_off
    end if

    ! nunk_tot = total amount of unknowns
    if ( present(prob_list_nunk) ) then
       nunk_tot=prob_list_nunk(1)
       do ielem=2, fspac%g_mesh%nelem
          nunk_tot=max(nunk_tot,prob_list_nunk(ielem))
       end do
    else
       nunk_tot = prob_nunk
    end if

    ! Loop over elements
    do ielem = 1, fspac%g_mesh%nelem
       ! Assign type of problem to ielem
       if ( present(prob_code) .and. present(prob_nunk) ) then
          fspac%lelem(ielem)%prob = prob_code
          nunk = prob_nunk
       else if ( present(prob_list) .and. present(prob_list_nunk) ) then
          fspac%lelem(ielem)%prob = prob_list(ielem)
          nunk = prob_list_nunk(ielem)
       else
          write(6,*) __FILE__,__LINE__,'ERROR:: problem not provided!!'
          assert(0==1)
       end if

       ! One only interpolation per element
       fspac%lelem(ielem)%nint = 1

       ! All the variables have the same interpolation
       call memalloc(nunk_tot,fspac%lelem(ielem)%iv,__FILE__,__LINE__)
       fspac%lelem(ielem)%iv(:) = 1

       ! Allocate the fixed info array
       call memalloc(1,fspac%lelem(ielem)%f_inf,__FILE__,__LINE__)

       ! Assign pointer to interpolation fixed information
       v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order(ielem)
       call fspac%ht_elem_info%put(key=v_key,val=fspac%cur_elinf,stat=istat)
       if ( istat == now_stored) then 
          ! Create fixed info if not constructed
          call fem_element_fixed_info_create(fspac%lelem_info(fspac%cur_elinf),f_type,              &
               &                             f_order(ielem),dim,created)
          assert(created)
          pos_elinf = fspac%cur_elinf
          fspac%cur_elinf = fspac%cur_elinf + 1
       else if ( istat == was_stored ) then
          call fspac%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          assert ( istat == key_found )
       end if
       fspac%lelem(ielem)%f_inf(1)%p => fspac%lelem_info(pos_elinf)

       ! Assign pointer to geometrical fixed information (assumed to be of order 1)
       v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       call fspac%ht_elem_info%put(key=v_key,val=fspac%cur_elinf,stat=istat)
       if ( istat == now_stored) then 
          ! Create fixed info if not constructed
          call fem_element_fixed_info_create(fspac%lelem_info(fspac%cur_elinf),f_type,              &
               &                             1,dim,created)
          assert(created)
          pos_elinf = fspac%cur_elinf
          fspac%cur_elinf = fspac%cur_elinf + 1
       else if ( istat == was_stored ) then
          call fspac%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          assert ( istat == key_found )
       end if
       fspac%lelem(ielem)%p_geo_info => fspac%lelem_info(pos_elinf)

       ! Compute number of DOFs in the elemental matrix associated to ielem
       if (fspac%kfl_scond == scond_on) then
          nnode = fspac%lelem(ielem)%f_inf(1)%p%nnode-fspac%lelem(ielem)%f_inf(1)%p%nodes_obj(dim+1)
       else
          nnode = fspac%lelem(ielem)%f_inf(1)%p%nnode
       end if
       fspac%lelem(ielem)%ldof = nunk*nnode

       ! Assign pointer to p_mat and p_vec in ielem
       call fspac%ht_pos_elmat%put(key=nunk*nnode,val=fspac%cur_elmat,stat=istat)
       call fspac%ht_pos_elvec%put(key=nunk*nnode,val=fspac%cur_elvec,stat=istat)
       if ( istat == now_stored ) then
          call elmat_alloc ( 1, 1, nunk*nnode,fspac%lelmat(fspac%cur_elmat), scal )
          pos_elmat = fspac%cur_elmat
          fspac%cur_elmat = fspac%cur_elmat + 1
          call elvec_alloc (1,nunk*nnode,fspac%lelvec(fspac%cur_elvec), scal )
          pos_elvec = fspac%cur_elvec
          fspac%cur_elvec = fspac%cur_elvec + 1
       else if ( istat == was_stored ) then
          call fspac%ht_pos_elmat%get(key=nunk*nnode,val=pos_elmat,stat=istat)
          assert ( istat == key_found )
          call fspac%ht_pos_elvec%get(key=nunk*nnode,val=pos_elvec,stat=istat)
          assert ( istat == key_found )
       end if       
       fspac%lelem(ielem)%p_mat => fspac%lelmat(pos_elmat)
       fspac%lelem(ielem)%p_vec => fspac%lelvec(pos_elvec)
       
       ! Allocate elem2dof, unkno, p_nod, integ and jvars
       call memalloc( nnode, nunk_tot, fspac%lelem(ielem)%elem2dof, __FILE__,__LINE__ )
       fspac%lelem(ielem)%elem2dof = 0
       call memalloc( nnode, nunk_tot, tdim,fspac%lelem(ielem)%unkno, __FILE__,__LINE__)
       fspac%lelem(ielem)%unkno = 0.0_rp
       call memalloc( nnode+1, fspac%lelem(ielem)%p_nod, __FILE__,__LINE__)
       fspac%lelem(ielem)%p_nod = 0
       call memalloc(1,fspac%lelem(ielem)%integ,__FILE__,__LINE__)
       call memalloc( nnode, nunk_tot, fspac%lelem(ielem)%jvars, __FILE__,__LINE__ )
       fspac%lelem(ielem)%jvars = 0

       ! Assign pointers to volume integration
       ltype(1) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order(ielem)
       ltype(2) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       v_key    = (max_ndime+1)*(max_FE_types+1)*(max_order) * ltype(1) + ltype(2)
       call fspac%ht_pos_vol_integ%put(key=v_key, val=fspac%cur_lvoli, stat = istat)
       if ( istat == now_stored ) then
          call integ_create(f_type,f_type,dim,1,f_order(ielem),fspac%lvoli(fspac%cur_lvoli),khie = khir)
          pos_voint       = fspac%cur_lvoli
          fspac%cur_lvoli = fspac%cur_lvoli + 1
       else if ( istat == was_stored ) then
          call fspac%ht_pos_vol_integ%get(key=v_key,val=pos_voint,stat=istat)
          assert ( istat == key_found )
       end if
       fspac%lelem(ielem)%integ(1)%p => fspac%lvoli(pos_voint)
       
    end do

  end subroutine fem_space_fe_list_create_one_int_p

  !==================================================================================================
  ! Fill the fem_space assuming that all elements are of type f_type but each variable has different
  ! interpolation order
  subroutine fem_space_fe_list_create_var_int( fspac, nint, iv, f_type, f_order, dim, tdim, &
       prob_code, prob_nunk, prob_list, prob_list_nunk, khie, scond  )
    implicit none
    type(fem_space), intent(inout),target  :: fspac
    integer(ip),     intent(in)       :: nint, iv(:), f_type, f_order(nint), dim, tdim
    integer(ip), optional, intent(in) :: prob_code, prob_nunk
    integer(ip), optional, intent(in) :: prob_list(:), prob_list_nunk(:)
    integer(ip), optional, intent(in) :: khie
    integer(ip), optional, intent(in) :: scond

    integer(ip) :: nunk, v_key, ltype(2), khir, nnode, max_nnode, nunk_tot
    integer(ip) :: ielem, istat, pos_elmat, pos_elinf,pos_elvec, pos_voint, iint, iunk, lndof
    logical(ip) :: created

    ! Hierarchical flag
    if (present(khie)) then
       khir = khie
    else
       khir = 0
    end if
    
    ! Static condensation flag
    if (present(scond)) then
       fspac%kfl_scond = scond
    else
       fspac%kfl_scond = scond_off
    end if

    ! nunk_tot = total amount of unknowns
    if ( present(prob_list_nunk) ) then
       nunk_tot=prob_list_nunk(1)
       do ielem=2, fspac%g_mesh%nelem
          nunk_tot=max(nunk_tot,prob_list_nunk(ielem))
       end do
    else
       nunk_tot = prob_nunk
    end if

    ! Loop over elements
    do ielem = 1, fspac%g_mesh%nelem
       ! Assign type of problem to ielem
       if ( present(prob_code) .and. present(prob_nunk) ) then
          fspac%lelem(ielem)%prob = prob_code
          nunk = prob_nunk
       else if ( present(prob_list) .and. present(prob_list_nunk) ) then
          fspac%lelem(ielem)%prob = prob_list(ielem)
          nunk = prob_list_nunk(ielem)
       else
          write(6,*) 'Error in fem_space_fe_list_create'
          assert(0==1)
       end if
       assert(size(iv,1)==nunk_tot)

       ! Set #interpolations per element
       fspac%lelem(ielem)%nint = nint

       ! Set interpolations x variable in iv
       call memalloc(nunk_tot,fspac%lelem(ielem)%iv,__FILE__,__LINE__)
       fspac%lelem(ielem)%iv = iv

       ! Allocate the fixed info array
       call memalloc(nint,fspac%lelem(ielem)%f_inf,__FILE__,__LINE__)

       ! Assign pointer to interpolation fixed information
       do iint=1,nint
          v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order(iint)
          call fspac%ht_elem_info%put(key=v_key,val=fspac%cur_elinf,stat=istat)
          if ( istat == now_stored) then 
             call fem_element_fixed_info_create(fspac%lelem_info(fspac%cur_elinf),f_type,              &
                  &                             f_order(iint),dim,created)
             assert(created)
             pos_elinf = fspac%cur_elinf
             fspac%cur_elinf = fspac%cur_elinf + 1
          else if ( istat == was_stored ) then
             call fspac%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
             assert ( istat == key_found )
          end if
          fspac%lelem(ielem)%f_inf(iint)%p => fspac%lelem_info(pos_elinf)
        
       end do

       ! Assign pointer to geometrical fixed information (assumed to be of order 1)
       v_key = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       call fspac%ht_elem_info%put(key=v_key,val=fspac%cur_elinf,stat=istat)
       if ( istat == now_stored) then 
          call fem_element_fixed_info_create(fspac%lelem_info(fspac%cur_elinf),f_type,              &
               &                             1,dim,created)
          assert(created)
          pos_elinf = fspac%cur_elinf
          fspac%cur_elinf = fspac%cur_elinf + 1
       else if ( istat == was_stored ) then
          call fspac%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
          assert ( istat == key_found )       
       end if
       fspac%lelem(ielem)%p_geo_info => fspac%lelem_info(pos_elinf)

       ! Compute number of DOFs in the elemental matrix associated to ielem
       lndof = 0
       max_nnode = 0
       do iunk=1,nunk_tot
          if (fspac%kfl_scond == scond_on) then
             nnode = fspac%lelem(ielem)%f_inf(iv(iunk))%p%nnode -                                   &
                  &  fspac%lelem(ielem)%f_inf(1)%p%nodes_obj(dim+1)
          else
             nnode = fspac%lelem(ielem)%f_inf(iv(iunk))%p%nnode 
          end if
          lndof = lndof + nnode
          max_nnode = max(max_nnode,nnode)
       end do
       fspac%lelem(ielem)%ldof = lndof

       ! Assign pointer to p_mat and p_vec in ielem
       call fspac%ht_pos_elmat%put(key=lndof,val=fspac%cur_elmat,stat=istat)
       call fspac%ht_pos_elvec%put(key=lndof,val=fspac%cur_elvec,stat=istat)
       if ( istat == now_stored ) then
          call elmat_alloc ( 1, 1, lndof,fspac%lelmat(fspac%cur_elmat), scal )
          pos_elmat = fspac%cur_elmat
          fspac%cur_elmat = fspac%cur_elmat + 1
          call elvec_alloc (1,lndof,fspac%lelvec(fspac%cur_elvec), scal )
          pos_elvec = fspac%cur_elvec
          fspac%cur_elvec = fspac%cur_elvec + 1
       else if ( istat == was_stored ) then
          call fspac%ht_pos_elmat%get(key=lndof,val=pos_elmat,stat=istat)
          assert ( istat == key_found )
          call fspac%ht_pos_elvec%get(key=lndof,val=pos_elvec,stat=istat)
          assert ( istat == key_found )
       end if
       fspac%lelem(ielem)%p_mat => fspac%lelmat(pos_elmat)
       fspac%lelem(ielem)%p_vec => fspac%lelvec(pos_elvec)
       
       ! Allocate elem2dof, unkno, p_nod, integ and jvars
       call memalloc( max_nnode, nunk_tot, fspac%lelem(ielem)%elem2dof, __FILE__,__LINE__ )
       fspac%lelem(ielem)%elem2dof = 0
       call memalloc( max_nnode, nunk_tot, tdim, fspac%lelem(ielem)%unkno, __FILE__,__LINE__)
       fspac%lelem(ielem)%unkno = 0.0_rp
       call memalloc( max_nnode+1, fspac%lelem(ielem)%p_nod, __FILE__,__LINE__)
       fspac%lelem(ielem)%p_nod = 0
       call memalloc(nint,fspac%lelem(ielem)%integ,__FILE__,__LINE__)
       call memalloc( max_nnode, nunk_tot, fspac%lelem(ielem)%jvars, __FILE__,__LINE__ )
       fspac%lelem(ielem)%jvars = 0

       ! Assign pointers to volume integration
       ltype(2) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)
       do iint = 1,nint
          ltype(1) = dim + (max_ndime+1)*f_type + (max_ndime+1)*(max_FE_types+1)*f_order(iint)
          v_key    = (max_ndime+1)*(max_FE_types+1)*(max_order) * ltype(1) + ltype(2)
          call fspac%ht_pos_vol_integ%put(key=v_key, val=fspac%cur_lvoli, stat = istat)
          if ( istat == now_stored ) then
             call integ_create(f_type,f_type,dim,1,f_order(iint),fspac%lvoli(fspac%cur_lvoli),     &
                  &            khie = khir,mnode=max_nnode)
             pos_voint       = fspac%cur_lvoli
             fspac%cur_lvoli = fspac%cur_lvoli + 1
          else if ( istat == was_stored ) then
             call fspac%ht_pos_vol_integ%get(key=v_key,val=pos_voint,stat=istat)
             assert ( istat == key_found )
          end if
          fspac%lelem(ielem)%integ(iint)%p => fspac%lvoli(pos_voint)

       end do
    end do

  end subroutine fem_space_fe_list_create_var_int

 !==================================================================================================
  subroutine dof_to_elnod( el2dof, idof, ivar, nobje, pnodob, lnode )
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: el2dof(:,:)
    integer(ip), intent(in)  :: idof, ivar, nobje
    integer(ip), intent(in)  :: pnodob(nobje+1)
    integer(ip), intent(out) :: lnode

    ! Locals
    integer(ip) :: iobje, jnode

    do iobje = 1,nobje
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

  subroutine fem_element_print ( lunou, elm )
    implicit none
    integer(ip)      , intent(in) :: lunou
    type(fem_element), intent(in) :: elm

    integer(ip) :: iint

    write (lunou, '(a)')     '*** begin finite element data structure ***'
    write (lunou, '(a,i10)') 'Problem: ', elm%prob    
    write (lunou,*) 'Element dofs: ', elm%elem2dof
    
    write (lunou,*) 'Number of interpolations: ', elm%nint
    write (lunou,*) 'Interpolation of each variable: ', elm%iv

    write (lunou,*) 'Fixed info of each interpolation: '
    do iint=1,elm%nint
       write (lunou, *) 'Type: ', elm%f_inf(iint)%p%ftype
       write (lunou, *) 'Order: ', elm%f_inf(iint)%p%order
       write (lunou, *) 'Nobje: ', elm%f_inf(iint)%p%nobje
       write (lunou, *) 'Nnode: ', elm%f_inf(iint)%p%nnode
       write (lunou, *) 'Nobje_dim: ', elm%f_inf(iint)%p%nobje_dim
       write (lunou, *) 'Nodes_obj: ', elm%f_inf(iint)%p%nodes_obj
       write (lunou, *) 'ndxob_i:  ', elm%f_inf(iint)%p%ndxob_i
       write (lunou, *) 'ndxob_j:  ', elm%f_inf(iint)%p%ndxob_j
       write (lunou, *) 'ntxob_i:  ', elm%f_inf(iint)%p%ntxob_i
       write (lunou, *) 'ntxob_j:  ', elm%f_inf(iint)%p%ntxob_j
       write (lunou, *) 'crxob_i:  ', elm%f_inf(iint)%p%crxob_i
       write (lunou, *) 'crxob_j:  ', elm%f_inf(iint)%p%crxob_j
!!$       write (lunou, *) 'Pointer object -> node: ', elm%f_inf(iint)%p%p_nodob
    end do

  end subroutine fem_element_print
  
  !==================================================================================================
  subroutine fem_space_free ( f )
    implicit none
    type(fem_space), intent(inout) :: f
    integer(ip)                    :: i,j
    
! lface_free
    if (allocated(f%lface)) then
       do i = 1, f%o_mesh%tface
          do j=1,f%lelem(f%lface(i)%nei_elem(1))%nint
             call array_free(f%lface(i)%o2n(j))
          end do
          call memfree(f%lface(i)%o2n,__FILE__,__LINE__)
       end do
       do i = 1, f%o_mesh%tface + f%o_mesh%bou_tface
          f%lface(i)%nei_elem  = 0
          f%lface(i)%pos_elem  = 0
          f%lface(i)%subface   = 0
          f%lface(i)%refinement_level = 0
          nullify (f%lface(i)%p_mat)
          nullify (f%lface(i)%p_vec)
          if (allocated(f%lface(i)%integ)) call memfree(f%lface(i)%integ,__FILE__,__LINE__)
          !nullify (f%lface(i)%aux_mat(1)%p)
          !nullify (f%lface(i)%aux_mat(2)%p)
       end do
      
       deallocate(f%lface)
       ! Face integrators free
       do i = 1,f%cur_lfaci-1 
          call integ_free( f%lfaci(i) )
       end do
       deallocate ( f%lfaci )
       f%cur_lfaci = 0
       call f%ht_pos_face_integ%free
    end if

    do i = 1, f%g_mesh%nelem
       f%lelem(i)%prob = 0
       call memfree( f%lelem(i)%elem2dof,__FILE__,__LINE__)
       call memfree( f%lelem(i)%unkno,__FILE__,__LINE__)
       call memfree( f%lelem(i)%jvars,__FILE__,__LINE__)
       if (allocated(f%lelem(i)%bc_value)) then
          do j=1,size(f%lelem(i)%bc_value) 
             if (allocated(f%lelem(i)%bc_value(j)%a)) call array_free(f%lelem(i)%bc_value(j))
          end do
          call memfree(f%lelem(i)%bc_value,__FILE__,__LINE__)
       end if
       nullify ( f%lelem(i)%p_geo_info )
       nullify ( f%lelem(i)%p_mat )
       nullify ( f%lelem(i)%p_vec )
       if(allocated(f%lelem(i)%f_inf)) call memfree(f%lelem(i)%f_inf,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%integ)) call memfree(f%lelem(i)%integ,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%iv))    call memfree(f%lelem(i)%iv   ,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%p_nod)) call memfree(f%lelem(i)%p_nod,__FILE__,__LINE__)
    end do
    deallocate ( f%lelem )

    do i = 1,f%cur_elmat-1
       call elmat_free( f%lelmat(i) )
    end do
    deallocate ( f%lelmat )

    do i = 1,f%cur_elvec-1
       call elvec_free( f%lelvec(i) )
    end do
    deallocate ( f%lelvec )

    f%cur_elmat = 0
    f%cur_elvec = 0

    call f%ht_pos_elmat%free
    call f%ht_pos_elvec%free

    do i = 1,f%cur_lvoli-1
       call integ_free( f%lvoli(i) )
    end do
    deallocate ( f%lvoli )
    f%cur_lvoli = 0
    call f%ht_pos_vol_integ%free

    do i = 1,f%cur_elinf-1
       call fem_element_fixed_info_free (f%lelem_info(i))
    end do
    deallocate (f%lelem_info)
    f%cur_elinf = 0
    call f%ht_elem_info%free

    nullify ( f%o_mesh )
    nullify ( f%g_mesh )

  end subroutine fem_space_free

  subroutine get_p_faces ( femsp, mesh )
    implicit none
    type(fem_space), intent(in)  :: femsp
    type(fem_mesh), intent(inout)    :: mesh
    
    integer(ip) :: i
    
    call memalloc ( mesh%nelem, mesh%p_face, __FILE__,__LINE__ )
    
    do i=1,mesh%nelem
          mesh%p_face(i) = mesh%pnods(i) + femsp%lelem(i)%f_inf(1)%p%nobje_dim(mesh%ndime) -1
    end do
    
  end subroutine get_p_faces

end module fem_space_names
