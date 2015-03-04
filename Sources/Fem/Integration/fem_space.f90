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
module fem_space_class
  ! Modules
  use types
  use memor
  use array_class
  use fem_triangulation_class
  use hash_table_class
  use integration_class
  use fem_space_types
#ifdef memcheck
  use iso_c_binding
#endif

  implicit none
# include "debug.i90"

  private

  integer(ip), parameter       :: max_global_interpolations  = 50   ! Maximum number of interpolations

  ! Information of each element of the FE space
  type fem_element

     integer(ip)                   :: prob          ! Problem to be solved
     integer(ip)                   :: nint          ! Number of different interpolations per element
     integer(ip)                   :: ldof          ! Dofs per element
     integer(ip),      allocatable :: material(:)   ! Material ! SB.alert : material can be used as p 
     logical(lg),      allocatable :: continuity(:)


     integer(ip)     , allocatable   :: elem2dof(:,:)   ! Map from elem to dof
     integer(ip)     , allocatable   :: iv(:)           ! Interpolation of each variable
     type(list_pointer), allocatable :: nodes_object(:) ! Nodes per object (including interior) (inter,obj)

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
     integer(ip)                        :: num_materials    ! Number of materials (maximum value)
     integer(ip)                        :: kfl_scond        ! Flag for static condensation 
     !(0=deactivated,1=activated)
     ! Elements and faces
     type(fem_element)    , allocatable :: lelem(:)         ! List of FEs
     type(fem_face)       , allocatable :: lface(:)         ! List of active faces

     type(fem_triangulation)  , pointer :: g_trian => NULL() ! Triangulation

     ! Elemental matrices and vector
     type (hash_table_ip_ip)            :: ht_pos_elmat
     type (hash_table_ip_ip)            :: ht_pos_elvec
     type(array_rp2)                    :: lelmat(max_global_interpolations)
     type(array_rp1)                    :: lelvec(max_global_interpolations)
     type(list)                         :: l_nodes_object(max_global_interpolations)
     integer(ip)                        :: cur_elmat
     integer(ip)                        :: cur_elvec

     ! Integrator
     type (hash_table_ip_ip)            :: ht_pos_vol_integ
     type (hash_table_ip_ip)            :: ht_pos_face_integ
     type(vol_integ)                    :: lvoli(max_global_interpolations)        
     type(face_integ)                   :: lfaci(max_global_interpolations)
     integer(ip)                        :: cur_lvoli
     integer(ip)                        :: cur_lfaci

     ! Element common information
     type (hash_table_ip_ip)            :: ht_elem_info
     type (fem_fixed_info), allocatable :: lelem_info(:)
     integer(ip)                         :: cur_elinf

  end type fem_space

  ! Memalloc interfaces
  interface memalloc
     module procedure memalloc_fem_fixed_info_p1, memalloc_fem_fixed_info_p2,                       &
          &           memalloc_fem_fixed_info_p3, memalloc_fem_fixed_info_p4,                       &
          &           memalloc_vol_integ_p1, memalloc_vol_integ_p2,                                 &
          &           memalloc_vol_integ_p3, memalloc_vol_integ_p4,                                 &
          &           memalloc_face_integ_p1, memalloc_face_integ_p2,                               &
          &           memalloc_face_integ_p3, memalloc_face_integ_p4
  end interface memalloc
  interface memrealloc
     module procedure memrealloc_fem_fixed_info_p1, memrealloc_fem_fixed_info_p2,                   &
          &           memrealloc_fem_fixed_info_p3, memrealloc_fem_fixed_info_p4,                   &
          &           memrealloc_vol_integ_p1, memrealloc_vol_integ_p2,                             &
          &           memrealloc_vol_integ_p3, memrealloc_vol_integ_p4,                              &
          &           memrealloc_face_integ_p1, memrealloc_face_integ_p2,                             &
          &           memrealloc_face_integ_p3, memrealloc_face_integ_p4
  end interface memrealloc
  interface memfree
     module procedure memfree_fem_fixed_info_p1, memfree_fem_fixed_info_p2,                         &
          &           memfree_fem_fixed_info_p3, memfree_fem_fixed_info_p4,                         &
          &           memfree_vol_integ_p1, memfree_vol_integ_p2,                                   &
          &           memfree_vol_integ_p3, memfree_vol_integ_p4,                                   &
          &           memfree_face_integ_p1, memfree_face_integ_p2,                                   &
          &           memfree_face_integ_p3, memfree_face_integ_p4
  end interface memfree
  interface memmovealloc
     module procedure memmovealloc_fem_fixed_info_p1, memmovealloc_fem_fixed_info_p2,               &
          &           memmovealloc_fem_fixed_info_p3, memmovealloc_fem_fixed_info_p4,               &
          &           memmovealloc_vol_integ_p1, memmovealloc_vol_integ_p2,                         &
          &           memmovealloc_vol_integ_p3, memmovealloc_vol_integ_p4,                         &
          &           memmovealloc_face_integ_p1, memmovealloc_face_integ_p2,                         &
          &           memmovealloc_face_integ_p3, memmovealloc_face_integ_p4
  end interface memmovealloc

  ! Types
  public :: fem_space, fem_element, fem_face,                      &
       &     memalloc, memrealloc, memfreep, memmovealloc

  ! Functions
  public :: fem_space_create, fem_space_fe_list_create, fem_element_print, fem_space_free,          &
       get_p_faces

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

  !**************************************************************************************************
  !**************************************************************************************************
  ! Specialization to allocate type(fem_fixed_info_pointer)
  !**************************************************************************************************
  !**************************************************************************************************
# define var_type type(fem_fixed_info_pointer)
# define var_size 8
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_fem_fixed_info_p1    
# define generic_memalloc_2      memalloc_fem_fixed_info_p2    
# define generic_memalloc_3      memalloc_fem_fixed_info_p3    
# define generic_memalloc_4      memalloc_fem_fixed_info_p4    
# define generic_memrealloc_1    memrealloc_fem_fixed_info_p1  
# define generic_memrealloc_2    memrealloc_fem_fixed_info_p2  
# define generic_memrealloc_3    memrealloc_fem_fixed_info_p3  
# define generic_memrealloc_4    memrealloc_fem_fixed_info_p4  
# define generic_memfree_1       memfree_fem_fixed_info_p1     
# define generic_memfree_2       memfree_fem_fixed_info_p2     
# define generic_memfree_3       memfree_fem_fixed_info_p3     
# define generic_memfree_4       memfree_fem_fixed_info_p4     
# define generic_memmovealloc_1  memmovealloc_fem_fixed_info_p1
# define generic_memmovealloc_2  memmovealloc_fem_fixed_info_p2
# define generic_memmovealloc_3  memmovealloc_fem_fixed_info_p3
# define generic_memmovealloc_4  memmovealloc_fem_fixed_info_p4
# include "memor.i90"

  !**************************************************************************************************
  !**************************************************************************************************
  ! Specialization to allocate type(vol_integ_pointer)
  !**************************************************************************************************
  !**************************************************************************************************
# define var_type type(vol_integ_pointer)
# define var_size 8
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_vol_integ_p1    
# define generic_memalloc_2      memalloc_vol_integ_p2    
# define generic_memalloc_3      memalloc_vol_integ_p3    
# define generic_memalloc_4      memalloc_vol_integ_p4    
# define generic_memrealloc_1    memrealloc_vol_integ_p1  
# define generic_memrealloc_2    memrealloc_vol_integ_p2  
# define generic_memrealloc_3    memrealloc_vol_integ_p3  
# define generic_memrealloc_4    memrealloc_vol_integ_p4  
# define generic_memfree_1       memfree_vol_integ_p1     
# define generic_memfree_2       memfree_vol_integ_p2     
# define generic_memfree_3       memfree_vol_integ_p3     
# define generic_memfree_4       memfree_vol_integ_p4     
# define generic_memmovealloc_1  memmovealloc_vol_integ_p1
# define generic_memmovealloc_2  memmovealloc_vol_integ_p2
# define generic_memmovealloc_3  memmovealloc_vol_integ_p3
# define generic_memmovealloc_4  memmovealloc_vol_integ_p4
# include "memor.i90"


  !**************************************************************************************************
  !**************************************************************************************************
  ! Specialization to allocate type(face_integ_pointer)
  !**************************************************************************************************
  !**************************************************************************************************
# define var_type type(face_integ_pointer)
# define var_size 8
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_face_integ_p1    
# define generic_memalloc_2      memalloc_face_integ_p2    
# define generic_memalloc_3      memalloc_face_integ_p3    
# define generic_memalloc_4      memalloc_face_integ_p4    
# define generic_memrealloc_1    memrealloc_face_integ_p1  
# define generic_memrealloc_2    memrealloc_face_integ_p2  
# define generic_memrealloc_3    memrealloc_face_integ_p3  
# define generic_memrealloc_4    memrealloc_face_integ_p4  
# define generic_memfree_1       memfree_face_integ_p1     
# define generic_memfree_2       memfree_face_integ_p2     
# define generic_memfree_3       memfree_face_integ_p3     
# define generic_memfree_4       memfree_face_integ_p4     
# define generic_memmovealloc_1  memmovealloc_face_integ_p1
# define generic_memmovealloc_2  memmovealloc_face_integ_p2
# define generic_memmovealloc_3  memmovealloc_face_integ_p3
# define generic_memmovealloc_4  memmovealloc_face_integ_p4
# include "memor.i90"

  !==================================================================================================
  ! Allocation of variables in fem_space according to the values in g_trian
  subroutine fem_space_create( fspac, g_trian )
    implicit none
    type(fem_space)  , intent(inout)        :: fspac
    type(fem_triangulation)   , intent(in), target   :: g_trian           

    integer(ip) :: istat

    allocate(fspac%lelem( g_trian%num_elems), stat=istat)
    if(istat/=0) then 
       write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
       write(6,'(a,i20)') 'Error code: ', istat
       write(6,'(a)') 'Caller routine: fem_space::fem_space_create::lelem'
    end if

    !  Initialization of pointers to triangulation
    fspac%g_trian => g_trian

    ! Allocation of elemental matrix and vector parameters
    call fspac%ht_pos_elmat%init(ht_length)
    fspac%cur_elmat = 1
    !    allocate(fspac%lelmat(max_interpolations), stat=istat)
    if(istat/=0) then 
       write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
       write(6,'(a,i20)') 'Error code: ', istat
       write(6,'(a)') 'Caller routine: fem_space::fem_space_create::lelmat'
    end if
    call fspac%ht_pos_elvec%init(ht_length)
    fspac%cur_elvec = 1
    !    allocate(fspac%lelvec(max_interpolations), stat=istat)
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
    !    allocate(fspac%lvoli(max_interpolations), stat=istat)
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
  ! Fill the fem_space assuming that all elements are of type f_type but each variable has different
  ! interpolation order
  subroutine fem_space_fe_list_create( fspac, nint, iv, continuity, material, f_type, f_order, dim, tdim, &
       prob_code, prob_nunk, prob_list, prob_list_nunk, khie, scond  )
    implicit none
    type(fem_space), intent(inout),target  :: fspac
    integer(ip),     intent(in)       :: nint, iv(:), material(:), f_type, f_order(nint), dim, tdim
    integer(ip), optional, intent(in) :: prob_code, prob_nunk
    integer(ip), optional, intent(in) :: prob_list(:), prob_list_nunk(:)
    integer(ip), optional, intent(in) :: khie
    integer(ip), optional, intent(in) :: scond 
    logical(lg),     intent(in)       :: continuity(:)

    integer(ip) :: nunk, v_key, ltype(2), khir, nnode, max_nnode, nunk_tot
    integer(ip) :: ielem, istat, pos_elmat, pos_elinf, pos_elvec, pos_voint, iint, iunk, lndof
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
       do ielem=2, fspac%g_trian%num_elems
          nunk_tot=max(nunk_tot,prob_list_nunk(ielem))
       end do
    else
       nunk_tot = prob_nunk
    end if

    ! void nodes object
    fspac%l_nodes_object(1)%n = max_nobje
    call memalloc( max_nobje+1, fspac%l_nodes_object(1)%p, __FILE__, __LINE__ )
    fspac%l_nodes_object(1)%p = 1
    call memalloc( 0, fspac%l_nodes_object(1)%l, __FILE__, __LINE__ )


    ! Loop over elements
    do ielem = 1, fspac%g_trian%num_elems
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
       assert(size(continuity,1)==nunk_tot)
       assert(size(material,1)==nunk_tot)

       ! Set #interpolations per element
       fspac%lelem(ielem)%nint = nint

       ! Set interpolations x variable in iv
       call memalloc(nunk_tot,fspac%lelem(ielem)%iv,__FILE__,__LINE__)
       fspac%lelem(ielem)%iv = iv

       ! Set continuity per interpolation
       call memalloc(nunk_tot,fspac%lelem(ielem)%continuity,__FILE__,__LINE__)
       fspac%lelem(ielem)%continuity = continuity

       ! Set continuity per interpolation
       call memalloc(nunk_tot,fspac%lelem(ielem)%material,__FILE__,__LINE__)
       fspac%lelem(ielem)%material = material

       ! Allocate the fixed info array
       call memalloc(nint,fspac%lelem(ielem)%f_inf,__FILE__,__LINE__)

       ! Allocate nodes per object info
       allocate( fspac%lelem(ielem)%nodes_object(nint), stat=istat )
       check( istat==0 )

       ! Assign pointer to interpolation fixed information and nodes per object
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
          if ( continuity(iint) ) then
             fspac%lelem(ielem)%nodes_object(iint)%p => fspac%lelem_info(pos_elinf)%ndxob
          else 
             fspac%lelem(ielem)%nodes_object(iint)%p => fspac%l_nodes_object(1) ! SB.alert : Think about hdG
          end if
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
          call array_create ( lndof, lndof, fspac%lelmat(fspac%cur_elmat) )
          pos_elmat = fspac%cur_elmat
          fspac%cur_elmat = fspac%cur_elmat + 1
          call array_create ( lndof, fspac%lelvec(fspac%cur_elmat) )
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
       !call memalloc( max_nnode+1, fspac%lelem(ielem)%p_nod, __FILE__,__LINE__)
       !fspac%lelem(ielem)%p_nod = 0
       call memalloc(nint,fspac%lelem(ielem)%integ,__FILE__,__LINE__)
       !call memalloc( max_nnode, nunk_tot, fspac%lelem(ielem)%jvars, __FILE__,__LINE__ )
       !fspac%lelem(ielem)%jvars = 0

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

  end subroutine fem_space_fe_list_create

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
       write (lunou, *) 'ndxob%p:  ', elm%f_inf(iint)%p%ndxob%p
       write (lunou, *) 'ndxob%l:  ', elm%f_inf(iint)%p%ndxob%l
       write (lunou, *) 'ntxob%p:  ', elm%f_inf(iint)%p%ntxob%p
       write (lunou, *) 'ntxob%l:  ', elm%f_inf(iint)%p%ntxob%l
       write (lunou, *) 'crxob%p:  ', elm%f_inf(iint)%p%crxob%p
       write (lunou, *) 'crxob%l:  ', elm%f_inf(iint)%p%crxob%l
!!$       write (lunou, *) 'Pointer object -> node: ', elm%f_inf(iint)%p%p_nodob
    end do

  end subroutine fem_element_print

  !==================================================================================================
  subroutine fem_space_free ( f )
    implicit none
    type(fem_space), intent(inout) :: f
    integer(ip)                    :: i,j

    ! SB.alert : To be thought for the new structure triangulation
    ! lface_free
    ! if (allocated(f%lface)) then
    !    j = 1
    !    do i = 1, f%g_trian%num_objects
    !       if ( trian%objects(i)%dimension == f%g_trian%num_dims-1 ) then
    !          f%lface(i)%nei_elem  = 0
    !          f%lface(i)%pos_elem  = 0
    !          f%lface(i)%subface   = 0
    !          f%lface(i)%refinement_level = 0
    !          nullify (f%lface(i)%p_mat)
    !          nullify (f%lface(i)%p_vec)
    !          if (allocated(f%lface(i)%integ)) call memfree(f%lface(i)%integ,__FILE__,__LINE__)
    !          if ( trian%objects(i)%border /= -1) then
    !             do j=1,f%lelem(f%lface(i)%nei_elem(1))%nint
    !                call array_free(f%lface(i)%o2n(j))
    !             end do
    !             call memfree(f%lface(i)%o2n,__FILE__,__LINE__)
    !          end if
    !          j = j+1
    !       end if
    !    end do

    !    deallocate(f%lface)
    !    ! Face integrators free
    !    do i = 1,f%cur_lfaci-1 
    !       call integ_free( f%lfaci(i) )
    !    end do
    !    f%cur_lfaci = 0
    !    call f%ht_pos_face_integ%free
    ! end if

    do i = 1, f%g_trian%num_elems
       f%lelem(i)%prob = 0
       call memfree( f%lelem(i)%elem2dof,__FILE__,__LINE__)
       call memfree( f%lelem(i)%unkno,__FILE__,__LINE__)
       !call memfree( f%lelem(i)%jvars,__FILE__,__LINE__)
       if (allocated(f%lelem(i)%bc_value)) then
          do j=1,size(f%lelem(i)%bc_value) 
             if (allocated(f%lelem(i)%bc_value(j)%a)) call array_free(f%lelem(i)%bc_value(j))
          end do
          call memfree(f%lelem(i)%bc_value,__FILE__,__LINE__)
       end if
       nullify ( f%lelem(i)%p_geo_info )
       nullify ( f%lelem(i)%p_mat )
       nullify ( f%lelem(i)%p_vec )
       do j = 1, f%lelem(i)%nint
          nullify ( f%lelem(i)%nodes_object(j)%p )
       end do
       if(allocated(f%lelem(i)%f_inf)) call memfree(f%lelem(i)%f_inf,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%integ)) call memfree(f%lelem(i)%integ,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%iv))    call memfree(f%lelem(i)%iv   ,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%continuity))    call memfree(f%lelem(i)%iv   ,__FILE__,__LINE__)
       if(allocated(f%lelem(i)%material))    call memfree(f%lelem(i)%iv   ,__FILE__,__LINE__)
       !if(allocated(f%lelem(i)%p_nod)) call memfree(f%lelem(i)%p_nod,__FILE__,__LINE__)
    end do
    deallocate ( f%lelem )

    do i = 1,f%cur_elmat-1
       call array_free( f%lelmat(i) )
    end do
    !deallocate ( f%lelmat )

    do i = 1,f%cur_elvec-1
       call array_free( f%lelvec(i) )
    end do

    !deallocate ( f%lelvec )

    f%cur_elmat = 0
    f%cur_elvec = 0

    call f%ht_pos_elmat%free
    call f%ht_pos_elvec%free

    do i = 1,f%cur_lvoli-1
       call integ_free( f%lvoli(i) )
    end do
    !deallocate ( f%lvoli )
    f%cur_lvoli = 0
    call f%ht_pos_vol_integ%free

    do i = 1,f%cur_elinf-1
       call fem_element_fixed_info_free (f%lelem_info(i))
    end do
    deallocate (f%lelem_info)
    f%cur_elinf = 0
    call f%ht_elem_info%free

    nullify ( f%g_trian )

  end subroutine fem_space_free

  subroutine get_p_faces ( femsp, trian )
    implicit none
    type(fem_space), intent(in)  :: femsp
    type(fem_triangulation), intent(inout)    :: trian

    !   integer(ip) :: i


    !   num_faces

    !   call memalloc ( trian%num_elems, num_faces, __FILE__,__LINE__ )

    !   do i=1,mesh%nelem
    !         mesh%p_face(i) = mesh%pnods(i) + femsp%lelem(i)%f_inf(1)%p%nobje_dim(mesh%ndime) -1
    !   end do

  end subroutine get_p_faces

  ! SB.alert : to be thought now

end module fem_space_class
